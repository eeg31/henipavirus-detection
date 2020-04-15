# Author: Emma E. Glennon
# Contributers to data and analysis: Nancy John and Olivier Restif
# Last edited: April 2020
# Purpose: estimate number of undetected spillover events and small
#          outbreaks based on estimated offspring distributions for 
#          Nipah virus in people and Hendra virus in horses
# Contact: eeg31@cam.ac.uk/eeglennon@gmail.com
#
#*******************************************************************************

require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)
require(scales)

set.seed(06042020)

# Settings (which parts of the analysis to run)

from_scratch <- TRUE
fit_models <- TRUE 

# Additional control parameters and constants

nsim <- 10000     #number of simulations
max_gen <- 50       #maximum number of generations for which to simulate an outbreak
max_index <- 20


fit_nbinom <- function(data, R0=NULL){

  if(!is.null(R0)) {
    lik <- function(k) sum(-dnbinom(data, mu=R0, size=k, log=TRUE))
    return(optim(c(k=0.5), lik, method="Brent", lower=0, upper=10000)$par)
  }
  
  lik <- function(par) {
    k <- par[1]
    R0 <- par[2]
    sum(-dnbinom(data, mu=R0, size=k, log=TRUE))
  }
  
  optim(c(k=0.5, R0=1), lik)
}


#setup and data loading*********************************************************
if(from_scratch) fit_models <- TRUE

# LOAD DATA
nipah <- rio::import('data/henipa/henipa-outbreaks.csv') %>%
         as_tibble %>%
         filter(virus=='NiV-B') #exclude outbreaks outside South Asia

hendra <- rio::import('data/henipa/henipa-outbreaks.csv') %>%
          as_tibble %>%
          filter(virus=='HeV')

nipah_dist <- rio::import('data/henipa/nipah-transmissions.xlsx') %>%
              as_tibble %>%
              mutate(count_min=count,
                     count_max=count
                     )
  
hendra_ll <- rio::import('data/henipa/hendra-linelist.csv') %>%
             as_tibble %>%
             mutate(possible_vals=transmission_max-transmission_min+1)

hendra_dist <- tibble(virus='HeV',
                      infected=0:max(hendra_ll$transmission_max, na.rm=TRUE),
                      source='Hendra linelist') %>%
               mutate(count_min=purrr::map_dbl(infected, function(i){
                        nrow(filter(hendra_ll, transmission_max==i & transmission_min==i))
                      }),
                      count_max=purrr::map_dbl(infected, function(i){
                        nrow(filter(hendra_ll, transmission_max>=i & transmission_min<=i))
                      }),
                      count=purrr::map_dbl(infected, function(i){
                        d <- filter(hendra_ll, transmission_max>=i & transmission_min<=i)
                        sum(1/d$possible_vals)
                      })
               )

secondary_dists <- full_join(hendra_dist, nipah_dist) %>%
                   mutate(density=purrr::map2_dbl(count, virus, function(i, v){
                          i/sum(filter(., virus==v)$count)
                   }))

ggplot(secondary_dists)+
  geom_line(aes(x=infected, y=density, color=virus, group=virus))+
  scale_x_discrete('number of secondary infections')+
  theme_classic()+
  scale_color_viridis_d()+
  facet_grid(rows='virus')+
  guides(color=FALSE)

#****************************************************************************************
#estimate Nipah index cases


nipah_index <- tibble(size=1:max(nipah$est_indeces)) %>%
  mutate(x=purrr::map_dbl(size, function(i) sum(nipah$est_indeces==i)),
         x_dens=x/sum(x)) %>%
  filter(size <= max_index)

d <- rep(nipah_index$size, nipah_index$x) - 1

lik_surface <- as_tibble(expand_grid(mu=seq(1, 20, 0.01),
                                     k=seq(0,5,0.01))) %>%
  mutate(nloglik = purrr::map2_dbl(mu, k, function(m, i){
    sum(-dnbinom(d, mu=m, size=i, log=TRUE))
  }),
  lik = exp(-nloglik),
  weight = lik/max(lik),
  id=1:nrow(.)
  )

sample_dists <- tibble(id = sample(lik_surface$id, size=nsim, 
                                   prob=lik_surface$weight, replace=TRUE)) %>%
                mutate(sim_id=1:nrow(.)) %>%
  left_join(select(lik_surface, c('id','mu','k','weight'))) %>%
  expand_grid(size=1:max_index) %>%
  mutate(row=1:nrow(.),
         prob = purrr::map_dbl(row, function(i){
           dnbinom(size[i] - 1, mu=mu[i], size=k[i])
         })
  )
nipah_index <- sample_dists

g1 <- ggplot(lik_surface)+
  geom_raster(aes(x=mu, y=k, fill=weight))+
  theme_classic()+
  scale_x_continuous(expression(mu))+
  scale_y_continuous('k')+
  scale_fill_viridis_c('relative likelihood')

g2 <- ggplot(sample_dists) +
  geom_line(aes(x=size, y=prob, group=id), alpha=0.01) +
  theme_classic()+
  scale_x_continuous('number of index cases per cluster')+
  scale_y_continuous('density')

cowplot::plot_grid(g1, g2, nrow=2)

#****************************************************************************************
if(from_scratch){
  
  run_params <- tibble(virus=c('NiV', 'HeV')) %>%
                mutate(R0_mean=c(0.33,  # Nikolay 2019
                                 0.35), # estimated from linelist (median of R0_min and R0_max; see below)
                       R0_min=c(0.19,   # Nikolay 2019
                                0.15),  # estimated from linelist (mean of transmission_min)
                       R0_max=c(0.59,   # Nikolay 2019
                                0.5),   # estimated from linelist (mean of transmission_max)
                       row=1:nrow(.),
                       R0_sample=purrr::map(row, function(i){
                          range95 <- R0_max[i] - R0_min[i]
                          out <- rnorm(nsim, mean=R0_mean[i], sd=range95/3.92)
                          sapply(out, max, 0)
                       })) %>%
                unnest(cols='R0_sample') %>%
                mutate(k_sample=purrr::map2_dbl(R0_sample, virus, function(i, v){
                          x <- filter(secondary_dists, virus==v)
                          x <- sample(size=10000, x=x$infected, prob=x$density, replace=TRUE)
                          
                          fit_nbinom(data=x, R0=i)
                       }),
                       sim_id=1:nrow(.)
                       ) %>%
                left_join(nipah_index) %>%
                mutate(index=purrr::map2_dbl(mu, k, function(m, i) rnbinom(1, mu=m, size=i) + 1)
                )
  
  save(run_params, file='results/run_params.rda')
} else {
  load(file='results/run_params.rda')
}

#****************************************************************************************

# calculate and sample probability generating function for a stochastic branching process
# where random variable X~nBinom(R0, k)
sample_PGF <- function(R0, k, n_index=1){
  inner_PGF <- function(n_0) sum(rnbinom(n_0, mu=R0, size=k))
    
  gen <- 0
  out <- n_index
  n <- n_index
  while(gen < max_gen & n > 0){
    n <- inner_PGF(n)
    out <- out + n
    gen <- gen + 1
  }
    
  return(out)
}
  
sims <- run_params %>%
        mutate(size=purrr::map2_dbl(R0_sample, k_sample, function(i,j) {
          index <- sample(nipah_index$size, size=1, prob=nipah_index$est)
          sample_PGF(i, j, index)
        }))
  

#add observed outbreaks
ggplot(sims)+
  geom_point(aes(x=size, color=virus), stat='density')+
  scale_x_discrete('outbreak size')+
  facet_grid(rows='virus') +
  theme_classic()+
  scale_color_viridis_d()+
  guides(color=FALSE)


if(fit_models){
  #cluster apply to each row (R0 and k): fit observation function to outbreak size distribution
  
  
  fit_observation <- function(observed, expected, alpha=NULL, beta=NULL){
    expected <- expected[1:length(observed)]
    
    nlik <- function(par){
      alpha <- par['alpha']
      beta <- par['beta']
      
      f <- sapply(1:length(observed), function(i) (1+exp(beta-i))^(-alpha))
      
      sum(-dmultinom(observed, size=sum(observed), prob=f*expected, log=TRUE))
    }
    
    if(!is.null(alpha) & !is.null(beta)) return(nlik(par=c(alpha=alpha, beta=beta)))
    
    optim(par=c(alpha=5, beta=0.5), nlik)
  }
  
  NiV_counts <- tibble(size=1:max(nipah$cases)) %>%
                mutate(density=purrr::map_dbl(size, function(i){
                    (nrow(filter(sims, virus=='NiV' & size==i))+1)/nrow(filter(sims, virus=='NiV'))
                  }),
                  observed=purrr::map_dbl(size, function(i){
                    sum(nipah$cases==i)
                  })
                )
  
  lik_surface <- as_tibble(expand_grid(alpha=exp(seq(0,40,.1))-1,
                                       beta=seq(-40,10,0.1))) %>%
                 mutate(nloglik = purrr::map2_dbl(alpha, beta, function(a, b){
                   fit_observation(observed = NiV_counts$observed,
                                   expected = NiV_counts$density, 
                                   alpha=a, beta=b)
                   }),
                   lik = exp(-nloglik),
                   weight = lik/max(lik),
                   id=1:nrow(.)
                 )
  
  sample_dists <- tibble(id = sample(lik_surface$id, size=nsim, prob=lik_surface$weight, replace=TRUE)) %>%
                  left_join(select(lik_surface, c('id','alpha','beta','weight'))) %>%
                  expand_grid(size=1:20) %>%
                  mutate(row=1:nrow(.),
                         prob = purrr::map_dbl(row, function(i){
                            (1+exp(beta[i]-size[i]))^(-alpha[i])
                          })
                  )
  
  g1 <- ggplot(lik_surface)+
        geom_raster(aes(x=log(1+alpha), y=beta, fill=weight))+
        theme_classic()+
        scale_x_continuous(expression(alpha))+
        scale_y_continuous(expression(beta))+
        scale_fill_viridis_c('relative likelihood')
  
  g2 <- ggplot(sample_dists) +
        geom_line(aes(x=size, y=prob, group=id), alpha=0.1) +
        theme_classic()+
        scale_x_continuous('number of cases')+
        scale_y_continuous('observation probability')
        
  cowplot::plot_grid(g1, g2, nrow=2)
  
  
} else {
  #load data
}