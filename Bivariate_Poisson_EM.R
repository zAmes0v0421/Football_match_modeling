library(bivpois)
library(alabama)
library(dplyr)

# Probability of Bivariate Poisson = (x, y)
# dbp(x, y, lambda1 = lam1, lambda2 = lam2, lambda3 = lam3) lambda3 refers to the common

# The 1st is not used.
conditional_porb = function(k, x, y, lam) {
  dpois(x - k, lambda = lam[1]) *
    dpois(y - k, lambda = lam[2]) *
    dpois(k, lambda = lam[3]) /
    dbp(x, y, lambda1 = lam[1], lambda2 = lam[2], lambda3 = lam[3])
}

# Conditional probability 
s_e_step = function(x, y, lam1, lam2, lam3) {
  k = 0:min(x, y)
  num = dpois(x - k, lam1) * dpois(y - k, lam2) * dpois(k, lam3) # single value minus vector allowed in R
  den = sum(num)  # equals bivariate Poisson probability
  cond_prob = num / den
  sum(k * cond_prob)
}


# Data cleaning. Using 2025-2026 England Prime league (truncated to 12 teams)

matches = read.csv("D:/RStudio/File_Coding/Projects/E0.csv", stringsAsFactors = FALSE)

matches_s = matches[, 4:7]

teams_group = c("Nott'm Forest", "Crystal Palace", "Leeds", "Everton", "Newcastle", "Liverpool",
                "Fulham", "Chelsea", "Brentford", "Sunderland", "Brighton", "Bournemouth")

matches = matches_s %>% filter(HomeTeam %in% teams_group)
matches = matches %>% filter(AwayTeam %in% teams_group)

# Indicator (f and g mapping) matrix making
home_mat = matrix(0, nrow = 132, ncol = 12)
away_mat = matrix(0, nrow = 132, ncol = 12)
colnames(home_mat) = teams_group
colnames(away_mat) = teams_group

home_mat[cbind(1:132, match(matches$HomeTeam, teams_group))] = 1

for (i in 1:132) {
  away_mat[i, matches$AwayTeam[i]] = 1
}

# Score matrix making
scores_data = as.matrix(matches[, 3:4])



# The trick is to use R’s closures: the E-step function takes data and current parameters, 
# computes the expected sufficient statistics, and then returns a new function 
# whose only argument is the candidate parameter vector. 
# The data and expectations are captured inside that returned function, 
# so the M-step’s optimizer never sees them.

# E step
# iter_att_def and att_def contain mu(-1) and home(-2)
Q_maker = function(score_data, home_m, away_m, iter_att_def, iter_beta){
  n = nrow(score_data)
  p = length(iter_att_def) / 2 - 1
  iter_att = iter_att_def[1:p]
  iter_def = iter_att_def[(p+1):(2*p)]
  iter_home = iter_att_def[2*p+1]
  iter_mu = iter_att_def[2*p+2]
  
  s = c()
  x = score_data[, 1]
  y = score_data[, 2]
  
  # lambda calculation
  home_lambda = exp(home_m %*% iter_att + away_m %*% iter_def + iter_home + iter_mu)
  away_lambda = exp(home_m %*% iter_def + away_m %*% iter_att + iter_mu)
  co_lambda = exp(home_m %*% iter_beta)
  iter_lambda = cbind(home_lambda, away_lambda, co_lambda)
  
  # Compute for each s_i
  for (i in 1:n){
    s = c(s, s_e_step(score_data[i,1], score_data[i,2],
                      iter_lambda[i,1], iter_lambda[i,2], iter_lambda[i, 3]))
  }
  # Q function 
  Q_one = function(att_def){ # att_def contains mu(-1) and home(-2)
    att = att_def[1:p]
    def = att_def[(p+1):(2*p)]
    home = att_def[2*p+1]
    mu = att_def[2*p+2]
    
    lambda_one = exp(home_m %*% att + away_m %*% def + home + mu)
    lambda_two = exp(home_m %*% def + away_m %*% att + mu)
    sum((x - s) * log(lambda_one)) - sum(lambda_one) + sum((y - s) * log(lambda_two)) - sum(lambda_two)
  }
  
  Q_two = function(beta){
    lambda_three = exp(home_m %*% beta)
    sum(s * log(lambda_three)) - sum(lambda_three)
  }
  
  Q = list(Q_1 = Q_one, Q_2 = Q_two)
  return(Q)
}

# M step. 
Q_one_optimizer = function(Q, att_def_start){ # att_def_start contains mu(-1) and home(-2)
  p = length(att_def_start) / 2 - 1
  
  heq = function(att_def) {
    c(sum(att_def[1:p]),           # sum(lam1) = 0
      sum(att_def[(p+1):(2*p)]))   # sum(lam2) = 0
  }
  result = auglag(par = att_def_start, # auglag: minimization
                   fn = function(theta) -Q(theta), # For maximization
                   heq = heq)
  return(result$par)
} 

Q_two_optimizer = function(Q, beta_start) {
  opt = optim(
    par = beta_start,
    fn = function(b) -Q(b),
    method = "BFGS"
  )
  opt$par
}

# Run EM algorithm
critical_value = 1e-3
max_iter = 10000
# Start
home = 0
mu = log(mean(c(scores_data)))
att_def = rep(0, 24)
att_def = c(att_def, home, mu)
beta = rep(0, 12)

for (i in 1: max_iter){
  att_def_iter = att_def
  beta_iter = beta
  
  # E step
  Q = Q_maker(scores_data, home_mat, away_mat, att_def_iter, beta_iter)
  
  # M step
  att_def = Q_one_optimizer(Q$Q_1, att_def_iter)
  beta = Q_two_optimizer(Q$Q_2, beta_iter)
  
  # Convergence judgement
  conv = max(abs(att_def_iter - att_def), abs(beta_iter - beta))
  cat(
    "conv =", format(conv, scientific = TRUE),
    " critical =", format(critical_value, scientific = TRUE),
    "\n"
  )
  if (conv < critical_value){
    cat("Converged at iteration:", i, "\n")
    break
  }
}

att_def_final = round(att_def, 3)
beta_final = round(beta, 3)

# Display 
att = att_def_final[1:12]
def = att_def_final[13:24]
home = att_def_final[25]
mu = att_def_final[26]

results = data.frame(
  Teams = teams_group,
  attack = att, 
  defense = def,
  home_covariance = beta_final 
)

# Calculate the observed likelihood by estimates
# log lambda calculator
log_lambda = function(home_matrix, att, away_matrix, def, home, mu, beta){
  home_lambda = home_matrix %*% att + away_matrix %*% def + home + mu
  away_lambda = home_matrix %*% def + away_matrix %*% att + mu
  co_lambda = home_matrix %*% beta
  log_lambda = cbind(home_lambda, away_lambda, co_lambda)
}

match_log_lambda = log_lambda(home_mat, att, away_mat, def, home, mu, beta_final)
colnames(match_log_lambda) = c("log_lambda_1", "log_lambda_2", "log_lambda_common")

# log likelihood for all matches
log_likelihood = 0
for (i in 1:nrow(scores_data)){
  log_likelihood = log_likelihood + dbp(scores_data[i,1], scores_data[i,2],
                                        lambda = exp(match_log_lambda[i, ]), logged = TRUE)
}
AIC_BiPois = -2*log_likelihood + 2*(22+1+1+12) # Bivariate Poisson has AIC of 815.


