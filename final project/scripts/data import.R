#-------- packages --------
library(tidyverse)
library(spdep)
library(ggthemes)
library(ggmap)
library(viridis)
library(lubridate)
library(gstat)
library(sp)
library(sf)
library(classInt)
library(lmtest)
library(tseries)
library(broom)
library(mgcv)


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
				 month = month(pickup_datetime, label = TRUE))

taxi_data <- taxi_data %>%
	select(total_amount, pickup_longitude, pickup_latitude, passenger_count, trip_distance, wday, hour, month)

manhattan <- readRDS("manhattan.rds")

taxi_data_sample <- sample_n(taxi_data, 1000) %>%
	select(total_amount, pickup_longitude, pickup_latitude, passenger_count, trip_distance, wday, hour, month)

spdf <- SpatialPointsDataFrame(coords = taxi_data[, 2:3], data = taxi_data,
															 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

taxi_sf <- st_as_sf(spdf)

weights <- knn2nb(knearneigh(taxi_sf, k = 50))
weights <- nb2listw(weights, zero.policy = TRUE)

#-------- functions --------
linear_regression <- function(data, formula){
	lm(formula, data)
}

persp_color <- function(data, palette, n_colors = 100) {

	if(any(!is_numeric(data), na.rm = TRUE)) {
		stop("data contains non-numeric values.")
	}

	n_rows <- nrow(data)
	n_cols <- ncol(data)

	jet.colors <- colorRampPalette(palette)
	color <- jet.colors(n_colors)

	z = data

	zfacet <- z[-1, -1] + z[-1, -n_cols] + z[-n_rows, -1] + z[-n_rows, -n_cols]
	facetcol <- cut(zfacet, n_colors)

	return(color[facetcol])
}


#-------- visualization and summary--------
ggmap(manhattan, darken = 0.5) +
	scale_fill_viridis(option = 'plasma') +
	geom_bin2d(data = taxi_data_sample, aes(pickup_longitude, pickup_latitude), bins = 50, alpha = 0.6) +
	labs(x = "longitude",
			 y = "latitude",
			 fill = "count")

taxi_data %>%
	ggplot(aes(hour, total_amount, color = wday)) +
	geom_boxplot() +
	facet_wrap(~month) +
	theme_linedraw()


#-------- krigging --------
set.seed(1234)

# GAM
gam <- gam(total_amount ~ s(pickup_longitude) + s(pickup_latitude) + passenger_count + s(trip_distance) + wday + hour + month, family = "gaussian", data = taxi_data)

cbind(fitted(gam), residuals(gam)) %>%
	as_tibble() %>%
	ggplot(aes(V1, V2)) +
	geom_point()

# data
taxi_sf_complete <- bind_cols(taxi_sf, residuals(gam), fitted(gam)) %>%
	rename("residuals" = "...10",
				 "fitted" = "...11")

sample_size <- floor(0.80*nrow(taxi_sf_complete))

select <- sample(seq_len(nrow(taxi_sf_complete)), size = sample_size)

test <- taxi_sf_complete[select,]
validation <- taxi_sf_complete[-select,]

combined_data <- bind_rows(test %>%
													 	mutate(type = "test"),
													 validation %>%
													 	mutate(type = "validation"))

# variograms
taxi_variogram <- variogram(residuals ~ 1,
														data = select(taxi_sf_complete, -fitted),
														cutoff = 10)
plot(taxi_variogram)

spherical_fit <- fit.variogram(taxi_variogram, vgm(0.05, "Sph", 5, 0.05))
plot(taxi_variogram, spherical_fit)

# ordinary kriging
geosp <- as_Spatial(select(taxi_sf_complete, -fitted))
prediction_grid <- spsample(geosp, nrow(taxi_sf), type = "regular")

test_sp <- as_Spatial(select(test, -fitted))

predict <- as_Spatial(select(validation, -fitted))

spherical_krige <- krige(formula = residuals ~ 1, locations=geosp, model=spherical_fit,
												 newdata=prediction_grid, nmax=15)

plot_data <- spherical_krige %>% as_tibble() %>%
											 	mutate(fit = "spherical")

final_data <- bind_cols(head(taxi_sf_complete, -84), spherical_krige$var1.pred) %>%
	mutate(fitted = fitted + ...12) %>%
	select(-...12)

# validation
spherical_validation <- krige(total_amount ~ 1, locations=test_sp, newdata= validation, model=spherical_fit)
exponential_validation <- krige(total_amount ~ 1, locations=test_sp, newdata= validation, model=exponential_fit)

spherical_difference <- validation$total_amount - spherical_validation$var1.pred
exponential_difference <- validation$total_amount - exponential_validation$var1.pred

rmseSph <- sqrt(sum(spherical_difference^2)/length(spherical_difference))
MESph   <- sum(spherical_difference/length(spherical_difference))
rmseExp <- sqrt(sum(exponential_difference^2)/length(exponential_difference))
MEExp   <- sum(exponential_difference/length(exponential_difference))

# cross-validation
spherical_cross_validation <- krige.cv(ni ~ 1, test_sp, model=spherical_fit, nfold=nrow(test_sp))
exponential_cross_validation <- krige.cv(ni ~ 1, test_sp, model=exponential_fit, nfold=nrow(test_sp))

sphcvmean <- mean(spherical_cross_validation$residual)
expcvmean <- mean(exponential_cross_validation$residual)
sphcvrmse <- sqrt(mean(spherical_cross_validation$residual^2))
expcvrmse <- sqrt(mean(exponential_cross_validation$residual^2))
sphcvmsdr <- mean((spherical_cross_validation$residual^2)/ (spherical_cross_validation$var1.var))
expcvmsdr <- mean((exponential_cross_validation$residual^2)/ (exponential_cross_validation$var1.var))

ggplot(combined_data, aes(type, total_amount)) +
	geom_boxplot(fill = "gold", color = "black") +
	stat_summary(fun=mean, geom="point", shape = 21, size = 3, fill = "white") +
	theme_par() +
	theme(axis.title.x = element_blank())

ggplot(combined_data, aes(total_amount, color = type, fill = type)) +
	geom_histogram(alpha = 0.8, show.legend = F) +
	facet_wrap(~ type) +
	theme_par() +
	scale_color_ptol() +
	scale_fill_ptol()

tibble(SSE = c(attr(spherical_fit, "SSErr"), attr(exponential_fit, "SSErr")),
			 `variogram fit` = c("spherical", "exponential")) %>%
	knitr::kable(caption = "Sum of squared errors for 2 variogram model fits", digits = 5)

combined_data %>%
	filter(type == "test") %>%
	ggplot(aes(pickup_longitude, pickup_latitude)) +
	geom_bin2d(aes(color=total_amount), bins = 50, alpha = 0.6) +
	coord_equal() +
	scale_color_gradient_tableau() +
	theme_par() +
	labs(color = "total amount") +
	theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				legend.position = "bottom",
				legend.key.width=unit(1.5,"cm"),
				axis.text.y = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks = element_blank())

final_data %>%
	ggplot(aes(pickup_longitude, pickup_latitude)) +
	geom_point(aes(color=total_amount), show.legend = FALSE) +
	coord_equal() +
	scale_color_gradient_tableau() +
	theme_map()

plot_data %>%
	ggplot(aes(x1, x2)) +
	geom_tile(aes(fill=var1.var)) +
	coord_equal() +
	facet_wrap(~ fit) +
	scale_fill_gradient_tableau() +
	theme_par() +
	labs(fill = "variance") +
	theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				legend.position = "bottom",
				legend.key.width=unit(1.5,"cm"),
				axis.text.y = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks = element_blank())

summary(spherical_krige)$data %>%
	knitr::kable(caption = "Spherical kriging predicted values and variances", digits = 5)



#-------- autocorrelation --------
ggplot(final_data, aes(pickup_longitude, pickup_latitude, color = exp(fitted), alpha = exp(fitted)))+
	geom_point(size = 0.1, show.legend = F) +
	scale_color_gradient2_tableau() +
	coord_fixed(ratio = 1.1) +
	theme_map() +
	labs(fill = "fare") +
	theme(legend.position = "top")

ggplot(final_data) +
	geom_histogram(aes(exp(total_amount)), bins = 100, position = "dodge") +
	#geom_histogram(aes(exp(fitted)), bins = 100, alpha = 0.6, color = "red")
	theme_minimal()

model <- taxi_data %>% linear_regression(total_amount ~ pickup_longitude + pickup_latitude + hour + wday + month + trip_distance + trip_time_in_secs)

model %>%
	broom::tidy() %>%
	knitr::kable()

model %>%
	broom::glance() %>%
	knitr::kable()

bind_rows(bptest(model) %>% broom::tidy(),
					jarque.bera.test(model$residuals) %>% broom::tidy()) %>%
	knitr::kable()

final_data %>%
	mutate(residuals2 = total_amount - fitted) %>%
	ggplot() +
		geom_point(aes(x=fitted, y=residuals2)) +
		theme_minimal() +
		geom_hline(aes(yintercept = 0), color = "royalblue")

taxi_data %>%
	bind_cols(model$residuals) %>%
	rename("residuals" = "...21") %>%
	ggplot(aes(pickup_longitude, pickup_latitude, color = residuals, alpha = residuals))+
	geom_point(size = 0.1, show.legend = F) +
	viridis::scale_color_viridis() +
	coord_fixed(ratio = 1.1) +
	theme_map() +
	labs(fill = "residual") +
	theme(legend.position = "top")

moran.test(model$residuals, listw = weights, zero.policy = TRUE)  %>%
	broom::tidy() %>%
	knitr::kable()

moran <- moran.mc(model$residuals, nsim = 100, listw = weights, zero.policy = TRUE)

moran %>%
	broom::tidy() %>%
	knitr::kable()

moran %>%
	plot()

#-------- models --------
formula <- as.formula(total_amount ~ pickup_longitude + pickup_latitude + hour + wday + month + trip_distance + trip_time_in_secs)

model <- glm(formula, family = "gaussian", taxi_data)

lapply(lm.LMtests(model, listw = weights, zero.policy = TRUE, test = "all"), tidy) %>%
	bind_rows(.id = "test") -> lm_results

sp_error_model <- errorsarlm(model, data = taxi_sf, listw = weights, zero.policy = TRUE, na.action = na.omit)

sp_lag_model <- lagsarlm(formula, taxi_sf, listw = weights, zero.policy = TRUE)

car_model <- spautolm(model, listw = weights, family = "CAR", zero.policy = TRUE)

model_names <- c("OLS", "SARerr", "SARlag", "CAR")
models <- list(model, sp_error_model, sp_lag_model, car_model)
names(models) <- model_names

comparison <- tibble(model            = model_names,
										 `Log Likelihood` = sapply(models, logLik),
										 AIC              = sapply(models, AIC),
										 BIC              = sapply(models, BIC))

lapply(models, function(x) moran.test(residuals(x), weights, zero.policy = TRUE) %>% tidy()) %>%
	bind_rows(.id = "model") -> moran_results

lapply(models, function(x) jarque.bera.test(residuals(x)) %>% tidy()) %>%
	bind_rows(.id = "model") -> jarque_berra_results

lapply(models, function(x) shapiro.test(residuals(x)) %>% tidy()) %>%
	bind_rows(.id = "model") -> shapiro_results

bind_rows(models$OLS %>% bptest() %>% tidy(),
					models$SARerr %>% bptest.sarlm() %>% tidy(),
					models$SARlag %>% bptest.sarlm() %>% tidy()) %>%
	mutate(model = c("OLS", "SARerr", "SARlag")) %>%
	select(model, statistic, p.value, parameter) -> bp_results

model_outputs <- lapply(models, residuals) %>%
	bind_rows(.id = "model") %>%
	pivot_longer(cols = 2:191, names_to = "row", values_to = "residual") %>%
	left_join(lapply(models, fitted) %>%
							bind_rows(.id = "model") %>%
							pivot_longer(cols = 2:191, names_to = "row", values_to = "fitted"),
						by = c("model", "row"))

comparison %>% kable(caption = "Diagnostic result comparison", digits = 5)
lm_results[, -5] %>% kable(caption = "Lagrange multiplier diagnostics for spatial dependence", digits = 5)
jarque_berra_results[, -5] %>% kable(caption = "Jarque Bera test",
																		 digits = 5)

shapiro_results[, -4] %>% kable(caption = "Shapiro-Wilk normality test",
																digits = 5)

model_outputs %>%
	ggplot(aes(residual, color = model, fill = model)) +
	geom_histogram(position = "dodge", alpha = 0.7, show.legend = F) +
	facet_wrap(~model) +
	scale_color_ptol() +
	scale_fill_ptol() +
	theme_light()

bp_results %>% kable(caption = "Breusch-Pagan test", digits = 5)

model_outputs %>%
	ggplot(aes(fitted, residual, color = model, fill = model)) +
	geom_point(alpha = 0.7, show.legend = F) +
	facet_wrap(~model) +
	scale_color_ptol() +
	scale_fill_ptol() +
	theme_light()

moran_results[, -c(7:8)] %>% kable(caption = "Moran I test under randomisation", digits = 5)

model %>% tidy() %>% kable(caption = "Ordinary least squares model (OLS)", digits = 5)
sp_error_model %>% tidy() %>% kable(caption = "Simultaneous autoregressive error model (SARerr)", digits = 5)
sp_lag_model %>% tidy() %>% kable(caption = "Simultaneous autoregressive lag model (SARlag)", digits = 5)
summary(car_model)$Coef %>% kable(caption = "Conditionally autoregressive model (CAR)", digits = 5)


gls_model <- gls(total_amount ~ pickup_longitude + pickup_latitude, correlation = corSpher(form = ~pickup_longitude + pickup_latitude, nugget = TRUE), data = taxi_sf)

bind_rows(bptest(gls_model) %>% broom::tidy(),
					jarque.bera.test(gls_model$residuals) %>% broom::tidy()) %>%
	knitr::kable()

jarque.bera.test(gls_model$residuals)


soilm <- taxi_data_sample %>%
	select(row = pickup_longitude, column = pickup_latitude, total_amount) %>%
	pivot_wider(names_from = "column", values_from = "total_amount") %>%
	select(-row) %>%
	as.matrix()

trend <- medpolish(soilm)

fit <- matrix(nrow = length(x), ncol = length(y))

for (i in seq_along(x)){
	for (j in seq_along(y)){
		fit[i,j] <- trend$row[i] + trend$col[j] + trend$overall
	}}
