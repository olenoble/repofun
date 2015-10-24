require(plyr)
require(fArma)

#### Get Data
path_data = '~/Code/R/TrendRiskParity/data/'
allfiles = list.files(path_data)

output = data.frame()
for(f in allfiles){

	fut_names = strsplit(gsub(".csv","",f), "_")[[1]]
	fut_names = fut_names[!gsub("1", "xxx", fut_names) == fut_names][1]
	cat("Getting Data for ", fut_names, "\n", sep = "")

	temp = read.csv(paste(path_data, f, sep = ""), stringsAsFactors = FALSE)
	temp = subset(temp, !is.na(Settle))[,c("Date","Settle")]
	names(temp)[2] = fut_names
	temp$Date = as.Date(temp$Date)

	if (dim(output)[1] == 0){
		output = temp
	} else {
		output = merge(output, temp)
	}

}

output_ret = cbind(data.frame(Date = output$Date[-1]), apply(output[,-1], 2, function(x){-1 + exp(diff(log(x)))}))

#### Generate Schedule
day_of_month = 15
start_date = as.Date(paste(format(min(output$Date), "%Y-%m"), day_of_month, sep = "-"))
end_date = as.Date(paste(format(max(output$Date), "%Y-%m"), day_of_month, sep = "-"))

rolldates = seq(start_date, end_date, by = "1 month")
schedule = data.frame(roll = 0:(length(rolldates) - 2), start = rolldates[-length(rolldates)], end = rolldates[-1])
schedule = subset(schedule, end <= end_date)

monthly_return = ddply(schedule, .(roll), function(X,Y){return(apply(subset(Y, Date >= X$start & Date <= X$end)[,-1],2,sum))}, output_ret)
yearly_return = apply(monthly_return[,-1], 2, function(X,h){
  X = cumsum(X)
  N = length(X)
  X[h:N] - c(0,X[1:(N-h)])}, 12)
yearly_return = cbind(data.frame(roll = monthly_return$roll[12:dim(monthly_return)[1]]), yearly_return)

undl = "SF1"
plot(output$Date, output[[undl]], type = "l")

ret_m = monthly_return[,c("roll", undl)]; names(ret_m)[2] = "monthly_ret"
ret_m_prev = ret_m; ret_m_prev$roll = ret_m_prev$roll + 1; names(ret_m_prev)[2] = "prev_monthly_ret"
ret_m = merge(ret_m, ret_m_prev)

ret_y = yearly_return[,c("roll", undl)]; ret_y$roll = ret_y$roll + 1; names(ret_y)[2] = "yearly_ret"
ret_m = merge(ret_m, ret_y)
summary(lm(monthly_ret ~ prev_monthly_ret + yearly_ret, ret_m))


pengFit(output_ret[[undl]], method = "mean")@hurst$diag[2,]

#### returns & weighted returns
#### autocorrel
#### Hurst

###### Optim test
###
output_ret = cbind(data.frame(Date = output$Date[-1]), apply(output[,-1], 2, function(x){diff(log(x))}))

getCovar = function(X, start_date, end_date){
  X = subset(X, Date >= start_date & Date <= end_date)
  X$Date = NULL
  Y = apply(X, 2, function(x){x - mean(x)})
  return((t(Y) %*% Y) / dim(X)[1])
}


start_date = as.Date("2005-01-01")
end_date = as.Date("2005-07-01")

ratio = runif(dim(output_ret)[2]-1)*2-1
optim_ret = getCovar(output_ret, start_date, end_date)


obj_func = function(x, params){
  return(-sum(abs(params$ratio) * log(abs(x))))
}

obj_func_minvar = function(x, params){
  x = array(x, dim = c(1, length(x)))
  return(c(x %*% params$cov_mat %*% t(x)))
}

eq_const = function(x, params){
  return(sum(abs(x)))
}

ineq_const = function(x, params){
  x = array(x, dim = c(1, length(x)))
  return(sqrt(c(x %*% params$cov_mat %*% t(x))) - params$target_vol)
}

input_params = list(ratio = ratio, cov_mat = optim_ret, target_vol = 0.25)
x_init = rep(1, length(ratio)) / length(ratio) * sign(ratio)
x_lb = (sign(ratio) - 1) / 2
x_ub = (sign(ratio) + 1) / 2

## First check minimum variance
#require(nloptr)
# res_minvar = nloptr(x0 = rep(1, length(ratio)) / length(ratio), 
#          eval_f = obj_func_minvar,
#          lb = rep(-1, length(ratio)),
#          ub = rep(1, length(ratio)),
#          eval_g_eq = eq_const,
#          opts = list("algorithm" = "NLOPT_GN_ISRES", "xtol_rel" = 1.0e-8),
#          params = input_params
#   )

require('Rsolnp')
res_minvar = solnp(pars = x_init, 
        fun = obj_func_minvar,
        LB = x_lb,
        UB = x_ub,
        eqfun = eq_const,
        eqB = 1,
        params = input_params
        )

cat("Minimum volatility  = ", round(sqrt(obj_func_minvar(res_minvar$pars, input_params) * 252) * 100, 4), "%\n")

## Then run risk parity
input_params = list(ratio = ratio, cov_mat = optim_ret, target_vol = 1 / sqrt(252))
x_init = rep(1, length(ratio)) / length(ratio) * sign(ratio)
x_lb = (sign(ratio) - 1) / 2
x_ub = (sign(ratio) + 1) / 2

res_riskpar = solnp(pars = x_init,
       fun = obj_func,
       LB = x_lb,
       UB = x_ub,
       eqfun = eq_const,
       eqB = 1,
       ineqfun = obj_func_minvar,
       ineqLB = 0,
       ineqUB = input_params$target_vol ^ 2,
       params = input_params
      )


cat("Achieved volatility  = ", round(sqrt(obj_func_minvar(res_riskpar$pars, input_params) * 252) * 100, 4), "%\n")
