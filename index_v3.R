require(plyr)

load("output.dat")

output_ret = cbind(data.frame(Date = output$Date[-1]), apply(output[,-1], 2, function(x){diff(log(x))}))

#### Generate Schedule
day_of_month = 15
start_date = as.Date(paste(format(min(output$Date), "%Y-%m"), day_of_month, sep = "-"))
end_date = as.Date(paste(format(max(output$Date), "%Y-%m"), day_of_month, sep = "-"))

rolldates = seq(start_date, end_date, by = "1 month")
schedule = data.frame(roll = 0:(length(rolldates) - 2), start = rolldates[-length(rolldates)], end = rolldates[-1])
schedule = subset(schedule, end <= end_date)

monthly_return = ddply(schedule, .(roll), function(X,Y){return(apply(subset(Y, Date >= X$start & Date <= X$end)[,-1],2,sum))}, output_ret)


#### KF
fit_kf_params = function(X){
	kf_coef =  summary(lm(X[2:length(X)] ~ X[1:(length(X)-1)]))$coefficients[,1]
	return(list(mu = kf_coef[1], beta = kf_coef[2], last_ret = X[length(X)]))
}

generateKF_ts = function(X, params_input){
	N = length(X)
	P = params_input$P0
	out_kf = c()
	for(i in 1:N){
		m_test = params_input$coefs$mu
		b_test = params_input$coefs$beta
		H = array(c(1, params_input$coefs$last_ret), dim = c(1,2))
		
		x_prev = array(c(params_input$coefs$mu, params_input$coefs$beta), dim = c(2,1))
		P = P + params_input$Q		
		K = P %*% t(H) %*% solve(H%*%P%*%t(H) + params_input$R)
		
		out_kf = c(out_kf, H %*% x_prev)
		x_next = x_prev + K %*% (X[i] - H %*% x_prev)
		
		params_input$coefs$mu = x_next[1]
		params_input$coefs$beta = x_next[2]
		params_input$coefs$last_ret = X[i]
	}
	params_input$P = P
	params_input$sim = out_kf
	
	return(params_input)
}

params = list()
params$P0 = array(0, dim = c(2,2))
params$R = 0.01 / 1200
params$Q = diag(c(0.05, 0.05)/100)

params$coefs = fit_kf_params(monthly_return[0:11,"SPY"])
params_test = params
params_test$coefs$last_ret = monthly_return[1,"SPY"]

sim_test = generateKF_ts(monthly_return[2:11,"SPY"], params_test)

#####
require('Rsolnp')

params_calib = list()
params_calib$P0 = array(0, dim = c(2,2))
params_calib$ts_input = monthly_return[1:12,"SPY"]
params_calib$coefs = fit_kf_params(params_calib$ts_input)
params_calib$coefs$last_ret = params_calib$ts_input[1]

x_init = (c(0.05, 0.1, 0.1) / sqrt(12)) ** 2

obj_calib_kf = function(X, params_kf){
	params_kf$R = X[1]
	params_kf$Q = diag(X[2:3])
	
	res_sim = generateKF_ts(params_kf$ts_input, params_kf)
	err = sum((res_sim$sim - res_sim$ts_input) ** 2)
	return(err)
}


res_minvar = solnp(pars = x_init * 100, 
	fun = obj_calib_kf,
	LB = c(0,0,0),
    params_kf = params_calib
    )

#######
params = list()
params$P0 = array(0, dim = c(2,2))
params$R = res_minvar$pars[1] 
params$Q = diag(res_minvar$pars[2:3])
params$coefs = fit_kf_params(monthly_return[0:11,"SPY"])

params$coefs$beta = params$coefs$beta * -1

params$coefs$last_ret = monthly_return[1,"SPY"]

sim_test = generateKF_ts(monthly_return[,"SPY"], params)

ret_actual = cumprod(c(1, exp(monthly_return[,"SPY"]))) * 100
ret_sim = cumprod(c(1, exp(sim_test$sim))) * 100

plot(ret_actual)
lines(ret_sim, col = "red")

