#-------- packages --------
library(tidyverse)
library(ggthemes)
library(spdep)
library(ggmap)
library(viridis)
library(lubridate)
library(gstat)
library(sp)
library(sf)
library(lmtest)
library(tseries)
library(broom)
library(mgcv)
library(randomForest)

#-------- data and directory --------
paste0(here::here(), "/final project/data") %>% setwd()

fares <- read_csv("fares.csv")
trips <- read_csv("trips.csv")

taxi_data <- trips %>%
	left_join(fares, by = c("medallion", "hack_license", "vendor_id", "pickup_datetime")) %>%
	select(-hack_license, -vendor_id, -rate_code, -store_and_fwd_flag) %>%
	filter(total_amount > 0) %>%
	mutate(total_amount = log(total_amount)) %>%
	filter(pickup_latitude >= 40.70, pickup_latitude <= 40.83,
				 pickup_longitude >= -74.025, pickup_longitude <= -73.93) %>%
	mutate(hour = hour(pickup_datetime),
				 wday = wday(pickup_datetime, label = TRUE),
				 month = month(pickup_datetime, label = TRUE)) %>%
	select(total_amount, pickup_longitude, pickup_latitude, passenger_count, trip_distance, wday, hour, month)

taxi_data <- sample_n(taxi_data, 2000)

sample_rows <- sample(nrow(taxi_data), 0.75 * nrow(taxi_data))

train <- taxi_data[sample_rows, ]
test <- taxi_data[-sample_rows, ]

train_sf <- SpatialPointsDataFrame(coords = train[, 2:3], data = train,
															 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
test_sf <- SpatialPointsDataFrame(coords = test[, 2:3], data = test,
																	proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

train_sf <- st_as_sf(train_sf)
test_sf <- st_as_sf(test_sf)

manhattan <- readRDS("manhattan.rds")
#-------- visualization and summary--------
ggmap(manhattan, darken = 0.5) +
	scale_fill_viridis(option = 'plasma') +
	geom_bin2d(data = taxi_data, aes(pickup_longitude, pickup_latitude), bins = 50, alpha = 0.6) +
	labs(x = "longitude",
			 y = "latitude",
			 fill = "count")

taxi_data %>%
	ggplot(aes(hour, total_amount, color = wday)) +
	geom_boxplot() +
	facet_wrap(~month) +
	theme_linedraw()

bind_rows(test %>%
						mutate(type = "test"),
					train %>%
						mutate(type = "train")) %>%
	ggplot(aes(type, total_amount)) +
	geom_boxplot(fill = "gold", color = "black") +
	stat_summary(fun=mean, geom="point", shape = 21, size = 3, fill = "white") +
	theme_par() +
	theme(axis.title.x = element_blank())

#-------- GAM --------
gam <- gam(total_amount ~ s(pickup_longitude) + s(pickup_latitude) + passenger_count + s(trip_distance) + wday + hour + month, family = "gaussian", data = train)

cbind(fitted(gam), residuals(gam)) %>%
	as_tibble() %>%
	ggplot(aes(V1, V2)) +
	geom_point()

# data with residuals
train_sf$resid <- residuals(gam)

test_sf$pred <- predict(gam, newdata = test)

test_sf %>%
	mutate(resid = total_amount - pred) %>%
	summarize(rmse = sqrt(mean(resid^2)))


#-------- krigging --------
set.seed(1234)

# variograms
taxi_variogram <- variogram(resid ~ 1,
														data = select(train_sf, -total_amount),
														cutoff = 10)
plot(taxi_variogram)

spherical_fit <- fit.variogram(taxi_variogram, vgm(0.05, "Sph", 3, 0.04))
plot(taxi_variogram, spherical_fit)

# ordinary kriging
train_sp <- as_Spatial(select(train_sf, -total_amount))
#prediction_grid <- spsample(train_sp, nrow(train_sf), type = "regular")
test_sp <- as_Spatial(select(test_sf, -pred))

spherical_krige <- krige(formula = resid ~ 1,
												 locations = train_sp,
												 model = spherical_fit,
												 newdata = test_sp,
												 nmax=15)

final_data <- bind_cols(test_sf, spherical_krige$var1.pred, spherical_krige$var1.var) %>%
	rename("variance" = "...12") %>%
	mutate(pred = pred + ...11) %>%
	select(- ...11)

jarque.bera.test(residuals(model))
jarque.bera.test(residuals(gam))
final_data %>%
	mutate(resid = pred - total_amount) %>%
	as_tibble() %>%
	select(resid) %>%
	pull() %>%
	jarque.bera.test()

shapiro.test(residuals(model))
shapiro.test(residuals(gam))
final_data %>%
	mutate(resid = pred - total_amount) %>%
	as_tibble() %>%
	select(resid) %>%
  pull() %>%
	shapiro.test()



final_data %>%
	mutate(resid = pred - total_amount) %>%
	summarize(rmse = sqrt(mean(resid^2)))

final_data %>%
	ggplot(aes(pickup_longitude, pickup_latitude)) +
	geom_point(aes(color=exp(pred), alpha = pred), show.legend = TRUE) +
	coord_equal() +
	scale_color_gradient_tableau() +
	theme_par() +
	labs(fill = "total earnings",
			 alpha = "log(total earnings)",
			 subtitle = "Manhattan",
			 x = "longitude",
			 y = "latitude")

ggplot(final_data, aes(pickup_longitude, pickup_latitude, color = exp(total_amount), alpha = total_amount))+
	geom_point(size = 2, show.legend = T) +
	scale_color_gradient_tableau() +
	coord_equal() +
	theme_par() +
	labs(fill = "total earnings",
			 alpha = "log(total earnings)",
			 subtitle = "Manhattan",
			 x = "longitude",
			 y = "latitude")

ggplot(final_data, aes(pickup_longitude, pickup_latitude, color = exp(variance)))+
	geom_point(size = 2, show.legend = T) +
	scale_color_gradient_tableau() +
	coord_equal() +
	theme_map() +
	labs(fill = "fare") +
	theme(legend.position = "top")

ggplot(final_data) +
	geom_histogram(aes(exp(total_amount)), bins = 100, color = "royalblue", fill = "royalblue") +
	geom_histogram(aes(exp(pred)), bins = 100, alpha = 0.6, color = "red", fill = "red") +
	theme_minimal() +
	labs(y = "",
			 x = "total earnings")


#-------- randomForests --------
rf_model <- randomForest(total_amount ~ ., data = train, ntree = 500)

test$pred <-  predict(rf_model, newdata = test, type = "class")

ggplot(test, aes(pickup_longitude, pickup_latitude, color = exp(pred), alpha = pred))+
	geom_point(size = 2, show.legend = T) +
	scale_color_gradient_tableau() +
	coord_equal() +
	theme_par() +
	labs(fill = "total earnings",
			 alpha = "log(total earnings)",
			 subtitle = "Manhattan",
			 x = "longitude",
			 y = "latitude")

test %>%
	mutate(residual = pred - total_amount) %>%
	summarize(rmse = sqrt(mean(residual^2)))

test %>%
	ggplot(aes(hour, pred, color = wday)) +
	geom_boxplot() +
	facet_wrap(~month) +
	theme_linedraw()

final_data %>%
	ggplot(aes(hour, pred, color = wday)) +
	geom_boxplot() +
	facet_wrap(~month) +
	theme_linedraw()

#-------- models --------
weights_test <- knn2nb(knearneigh(test_sf, k = 3))
weights_test <- nb2listw(weights_test, zero.policy = TRUE)

weights_train <- knn2nb(knearneigh(train_sf, k = 3))
weights_train <- nb2listw(weights_train, zero.policy = TRUE)

formula <- as.formula(total_amount ~ pickup_longitude + pickup_latitude + passenger_count + trip_distance + wday + hour + month)

model <- glm(formula, family = "gaussian", train)

lapply(lm.LMtests(model, listw = weights_train, zero.policy = TRUE, test = "all"), tidy) %>%
	bind_rows(.id = "test") -> lm_results

sp_error_model <- errorsarlm(model, data = train_sf, listw = weights_train, zero.policy = TRUE, na.action = na.omit)

sp_lag_model <- lagsarlm(formula, train_sf, listw = weights_train, zero.policy = TRUE)

car_model <- spautolm(model, train_sf, listw = weights_train, family = "CAR", zero.policy = TRUE, na.action = na.omit)

model_names <- c("OLS", "SARerr", "SARlag", "CAR")
models <- list(model, sp_error_model, sp_lag_model, car_model)
names(models) <- model_names

comparison <- tibble(model            = model_names,
										 `Log Likelihood` = sapply(models, logLik),
										 AIC              = sapply(models, AIC),
										 BIC              = sapply(models, BIC))

test$pred <- predict(model, newdata = test)

test %>%
	mutate(resid = total_amount - pred) %>%
	summarize(rmse = sqrt(mean(resid^2)))

train$pred <- predict(sp_error_model)

train %>%
	mutate(resid = total_amount - pred) %>%
	summarize(rmse = sqrt(mean(resid^2)))

train$pred <- spatialreg::predict.sarlm(sp_lag_model, newdata = test, listw = weights_test, zero.policy = TRUE)

train %>%
	mutate(resid = total_amount - pred) %>%
	summarize(rmse = sqrt(mean(resid^2)))

