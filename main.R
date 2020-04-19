# Author: Emma E. Glennon
# Contributers to data and analysis: Nancy John and Olivier Restif
# Last edited: April 2020
# Purpose: estimate number of undetected spillover events and small
#          outbreaks based on estimated offspring distributions for 
#          Nipah virus in people and Hendra virus in horses
# Contact: eeg31@cam.ac.uk/eeglennon@gmail.com
#
#*******************************************************************************

# Load dependencies ************************************************************

require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)
require(lhs)
require(scales)

# Set seed for reproducibility *************************************************

set.seed(06042020)

# Settings (which parts of the analysis to run) ********************************

from_scratch <- TRUE
simulate_outbreaks <- TRUE
fit_models <- TRUE

# Additional control parameters and constants **********************************

nsim <- 10000    # number of Monte Carlo simulations
max_gen <- 1000   # maximum number of generations for which to simulate an outbreak
max_index <- 30   # maximum number of index cases in a single outbreak

if(from_scratch) {
  fit_models <- FALSE         # fix conflicting settings
  simulate_outbreaks <- FALSE # fix conflicting settings
}

# Fit a negative binomial distribution to observed frequencies, 
# with or without a fixed mean *************************************************

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


# Load data ********************************************************************

# all recorded bat-to-human Nipah virus clusters in South Asia up to April 2020
nipah <- rio::import('data/henipa/henipa-outbreaks.csv') %>%
         as_tibble %>%
         filter(virus=='NiV-B') #exclude outbreaks outside South Asia

# all recorded bat-to-horse Hendra virus clusters in Australia up to April 2020
hendra <- rio::import('data/henipa/henipa-outbreaks.csv') %>%
          as_tibble %>%
          filter(virus=='HeV')

# estimated secondary case distribution of Nipah virus (Nikolay et al. 2019)
nipah_dist <- rio::import('data/henipa/nipah-transmissions.xlsx') %>%
              as_tibble %>%
              mutate(count_min=count,
                     count_max=count
                     )

# estimated linelist of all recorded Hendra virus infections in horses
hendra_ll <- rio::import('data/henipa/hendra-linelist.csv') %>%
             as_tibble %>%
             mutate(possible_vals=transmission_max-transmission_min+1)

# Transform Hendra linelist into secondary infection distribution **************

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

# plot and save secondary infection distributions ***************************************

ggplot(secondary_dists)+
  geom_line(aes(x=infected, y=density, color=virus, group=virus))+
  scale_x_continuous('number of secondary infections')+
  theme_classic()+
  scale_color_viridis_d(option='A', begin=0.3, end=0.7)+
  facet_grid(rows='virus')+
  guides(color=FALSE)
ggsave(file.path('figures','secondary_infections.pdf'), width=4, height=5)

# estimate and sample distribution of henipa index cases per case cluster ***************

nipah_index <- tibble(size=1:max(nipah$est_indeces)) %>%
  mutate(x=purrr::map_dbl(size, function(i) sum(nipah$est_indeces==i)),
         x_dens=x/sum(x)) %>%
  filter(size <= max_index)

d <- rep(nipah_index$size, nipah_index$x) - 1

lik_surface <- as_tibble(expand_grid(mu=seq(1, 10, 0.001),
                                     k=seq(0.001,5,0.001))) %>%
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
nipah_index <- sample_dists %>%
               mutate(virus='NiV')
sample_dists <- sample_dists %>%
  group_by(size) %>%
  summarise(p2.5=quantile(prob, 0.025),
            p25=quantile(prob, 0.25),
            p75=quantile(prob, 0.75),
            p97.5=quantile(prob, 0.975),
            p50=median(prob),
            min=min(prob),
            max=max(prob))

# Plot and save spillovers per case distribution estimates (Nipah) ***********************

g1 <- ggplot(lik_surface)+
  geom_raster(aes(x=mu, y=k, fill=weight))+
  theme_classic()+
  scale_x_continuous(expression(mu))+
  scale_y_continuous('k')+
  scale_fill_viridis_c('relative likelihood', option='A')

g2 <- ggplot(sample_dists) +
  geom_segment(aes(x=size, xend=size, group=size, y=min, yend=max), 
               alpha=0.5, lwd=1) +
  geom_segment(aes(x=size, xend=size, group=size, y=p2.5, yend=p97.5), 
               alpha=0.5, lwd=2) +
  geom_segment(aes(x=size, xend=size, group=size, y=p75, yend=p25), 
               alpha=0.5, lwd=4) +
  geom_point(aes(x=size, y=p50), size=3) +
  theme_classic()+
  scale_x_continuous('number of index cases per cluster')+
  scale_y_continuous('density')

print(sample_dists[1,])

#cowplot::plot_grid(g1, g2, nrow=2)


# estimate and sample distribution of Hendra index cases per case cluster ***************

hendra_index <- tibble(size=1:max(hendra$est_indeces)) %>%
  mutate(x=purrr::map_dbl(size, function(i) sum(hendra$est_indeces==i)),
         x_dens=x/sum(x)) %>%
  filter(size <= max_index)

d <- rep(hendra_index$size, hendra_index$x) - 1

lik_surface <- as_tibble(expand_grid(mu=seq(1, 2, 0.001),
                                     k=seq(0.001,5,.001))) %>%
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
           a <- dnbinom(size[i] - 1, mu=mu[i], size=k[i])
         })
  )

hendra_index <- sample_dists %>%
                mutate(virus='HeV')
sample_dists <- sample_dists %>%
  group_by(size) %>%
  summarise(p2.5=quantile(prob, 0.025),
            p25=quantile(prob, 0.25),
            p75=quantile(prob, 0.75),
            p97.5=quantile(prob, 0.975),
            p50=median(prob),
            min=min(prob),
            max=max(prob))

# Plot and save spillovers per case distribution estimates (Hendra) **********************

g3 <- ggplot(lik_surface)+
  geom_raster(aes(x=mu, y=k, fill=weight))+
  theme_classic()+
  scale_x_continuous(expression(mu))+
  scale_y_continuous('k')+
  scale_fill_viridis_c('relative likelihood', option='A')

g4 <- ggplot(sample_dists) +
  geom_segment(aes(x=size, xend=size, group=size, y=min, yend=max), 
               alpha=0.5, lwd=1) +
  geom_segment(aes(x=size, xend=size, group=size, y=p2.5, yend=p97.5), 
               alpha=0.5, lwd=2) +
  geom_segment(aes(x=size, xend=size, group=size, y=p75, yend=p25), 
               alpha=0.5, lwd=4) +
  geom_point(aes(x=size, y=p50), size=2) +
  theme_classic()+
  scale_x_continuous('number of index cases per cluster')+
  scale_y_continuous('density') 

print(sample_dists[1,])

cowplot::plot_grid(g2, g4, nrow=2, labels=c("A",'B'))
ggsave(file.path('figures','spills-per-cluster.pdf'), width=4, height=5)

# sample (or load) distribution parameters for each Monte Carlo simulation **************

index_df <- rbind(nipah_index,
                  hendra_index)

if(from_scratch){
  
  run_params <- tibble(virus=c('NiV', 'HeV')) %>%
                mutate(R0_mean=c(0.33,  # Nikolay 2019
                                 0.35), # estimated from linelist (median of R0_min, R0_max)
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
                left_join(index_df) %>%  # sample the number of index cases
                mutate(index=purrr::map2_dbl(mu, k, function(m, i) rnbinom(1, mu=m, size=i) + 1)
                )
  
  save(run_params, file='results/run_params.rda')
} else {
  
  load(file='results/run_params.rda')
}

# calculate and sample probability generating function for a stochastic branching process
# where random variable X~nBinom(R0, k) *************************************************
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
  
# sample probability generating function for each Monte Carlo simulation ****************

if(simulate_outbreaks){
  
  sims <- run_params %>%
          mutate(row=1:nrow(.),
                 index=purrr::map_int(virus, function(v){
                   d <- filter(index_df, virus==v)
                   sample(d$size, size=1, prob=d$prob)
                 }),
                 size=purrr::map_dbl(row, function(i) {
                   sample_PGF(R0_sample[i], k_sample[i], index[i])
                 }))
  
  save(sims, file=file.path('results',
                            paste('run_params-', nsim, '.rda', sep="")))
} else {
  
  load(file.path('results', paste('run_params-', nsim, '.rda', sep="")))
}
  
# add observed outbreaks ****************************************************************
nipah_plot <- nipah %>%
                mutate(virus='NiV')
  
ggplot()+
  geom_histogram(aes(x=cases, y=..count../sum(..count..)), binwidth=2, data=nipah_plot, fill='black') +
  geom_histogram(aes(x=cases, y=..count../sum(..count..)), binwidth=2, data=hendra, fill='black') +
  geom_area(aes(x=size, y=..count../sum(..count..), fill=virus, color=virus), 
            stat='count', lwd=0.5, alpha=0.5, data=filter(sims, virus=='NiV')) +
  geom_area(aes(x=size, y=..count../sum(..count..), fill=virus, color=virus), 
            stat='count', lwd=0.5, alpha=0.5, data=filter(sims, virus=='HeV')) +
  scale_x_continuous('outbreak size') +
  scale_y_continuous('density')+
  facet_grid(rows='virus', scales='free_y') +
  theme_classic() +
  scale_fill_viridis_d(option='A', begin=0.3, end=0.7) +
  scale_color_viridis_d(option='A', begin=0.3, end=0.7) +
  guides(color=FALSE)
ggsave(file.path('figures','observed_sizes.pdf'), width=4, height=5)

if(fit_models){
  
  NiV_counts <- tibble(size=1:max(nipah$cases)) %>%
                mutate(density=purrr::map_dbl(size, function(i){
                    (nrow(filter(sims, virus=='NiV' & size==i))+1)/nrow(filter(sims, virus=='NiV'))
                  }),
                  observed=purrr::map_dbl(size, function(i){
                    sum(nipah$cases==i)
                  }))

  save(NiV_counts, file=file.path('results', 
                             paste('nipah_results-', nsim, '.rda', sep="")))
  
  
  HeV_counts <- tibble(size=1:max(hendra$cases)) %>%
    mutate(density=purrr::map_dbl(size, function(i){
      (nrow(filter(sims, virus=='HeV' & size==i))+1)/nrow(filter(sims, virus=='HeV'))
    }),
    observed=purrr::map_dbl(size, function(i){
      sum(hendra$cases==i)
    }))
  
  save(HeV_counts, file=file.path('results', 
                                  paste('hendra_results-', nsim, '.rda', sep="")))
} else {
  load(file.path('results', paste('nipah_results-', nsim, '.rda', sep="")))
  load(file.path('results', paste('hendra_results-', nsim, '.rda', sep="")))
}


fit_observation <- function(observed, expected, alpha=NULL, beta=NULL, tol=1e-300){
  
  expected <- expected[1:length(observed)]
  
  nlik <- function(par){
    alpha <- par['alpha']
    beta <- par['beta']
    
    f <- sapply(1:length(observed), function(i) (1+exp(beta-i))^(-alpha))
    if(all(f<tol)) return(Inf)
    if(all(f==1) & any(f*expected!=observed)) return(Inf)

    sum(-dmultinom(observed, size=sum(observed), prob=f*expected, log=TRUE))
  }
  
  if(!is.null(alpha) & !is.null(beta)) return(nlik(par=c(alpha=alpha, beta=beta)))
  
  optim(par=c(alpha=5, beta=0.5), nlik)
}

lik_surface_NiV <- as_tibble(expand_grid(alpha=seq(0,10,.01),
                                         beta=seq(-20,20,0.01))) %>%
                   mutate(nloglik = purrr::map2_dbl(alpha, beta, function(a, b){
                     fit_observation(observed = NiV_counts$observed,
                                     expected = NiV_counts$density, 
                                     alpha=a, beta=b)
                     }),
                     lik = exp(-nloglik),
                     weight = lik/max(lik),
                     id=1:nrow(.)
                   )

lhs <- lhs::randomLHS(n=nsim, k=2)
NiV_posterior <- tibble(alpha=qunif(p=lhs[,1], min=0, max=10),
                        beta=qunif(p=lhs[,2], min=-20, max=20)) %>%
  mutate(nloglik = purrr::map2_dbl(alpha, beta, function(a, b){
    fit_observation(observed = NiV_counts$observed,
                    expected = NiV_counts$density, 
                    alpha=a, beta=b)
  }),
  lik = exp(-nloglik),
  weight = lik/max(lik),
  id=1:nrow(.)
  ) 

# resample
NiV_posterior <- NiV_posterior[sample(NiV_posterior$id, size=nsim, 
                                      prob=NiV_posterior$weight, replace=TRUE), ] %>%
                 expand_grid(size=1:20) %>%
                 mutate(row=1:nrow(.),
                        prob = purrr::map_dbl(row, function(i){
                          (1+exp(beta[i]-size[i]))^(-alpha[i])
                        })
                 ) %>%
                 group_by(size) %>%
                 summarise(p2.5=quantile(prob, 0.025),
                           p25=quantile(prob, 0.25),
                           p75=quantile(prob, 0.75),
                           p97.5=quantile(prob, 0.975),
                           p50=median(prob),
                           min=min(prob),
                           max=max(prob))
  
sample_dists_NiV_plot <- tibble(id = sample(lik_surface_NiV$id, size=nsim, prob=lik_surface_NiV$weight, replace=TRUE)) %>%
                        left_join(select(lik_surface_NiV, c('id','alpha','beta','weight'))) %>%
                        expand_grid(size=1:20) %>%
                        mutate(row=1:nrow(.),
                               prob = purrr::map_dbl(row, function(i){
                                  (1+exp(beta[i]-size[i]))^(-alpha[i])
                                })
                        ) %>%
      group_by(size) %>%
      summarise(p2.5=quantile(prob, 0.025),
                p25=quantile(prob, 0.25),
                p75=quantile(prob, 0.75),
                p97.5=quantile(prob, 0.975),
                p50=median(prob),
                min=min(prob),
                max=max(prob))

  
g1 <- ggplot(lik_surface_NiV)+
      geom_raster(aes(x=alpha, y=beta, fill=weight))+
      theme_classic()+
      scale_x_continuous(expression(alpha))+
      scale_y_continuous(expression(beta))+
      scale_fill_viridis_c('relative likelihood', option='A') +
      ggtitle('Nipah virus')
  
g2 <- ggplot(sample_dists_NiV_plot) +
      geom_segment(aes(x=size, xend=size, group=size, y=min, yend=max), 
                   alpha=0.5, lwd=1) +
      geom_segment(aes(x=size, xend=size, group=size, y=p97.5, yend=p2.5), 
                   alpha=0.5, lwd=2) +
      geom_segment(aes(x=size, xend=size, group=size, y=p75, yend=p25), 
                   alpha=0.5, lwd=4) +
      geom_point(aes(x=size, y=p50), size=3) +      theme_classic()+
      scale_x_continuous('number of cases')+
      scale_y_continuous('observation probability')
        
cowplot::plot_grid(g1, g2, nrow=2)
print(NiV_posterior[1,])


lik_surface_HeV <- as_tibble(expand_grid(alpha=seq(0,3,.001),
                                         beta=seq(-25,100,0.1))) %>%
  mutate(nloglik = purrr::map2_dbl(alpha, beta, function(a, b){
    fit_observation(observed = HeV_counts$observed,
                    expected = HeV_counts$density, 
                    alpha=a, beta=b)
  }),
  lik = exp(-nloglik),
  weight = lik/max(lik),
  id=1:nrow(.)
  )

sample_dists_HeV_vis <- tibble(id = sample(lik_surface_HeV$id, size=nsim, prob=lik_surface_HeV$weight, replace=TRUE)) %>%
  left_join(select(lik_surface_HeV, c('id','alpha','beta','weight'))) %>%
  expand_grid(size=1:20) %>%
  mutate(row=1:nrow(.),
         prob = purrr::map_dbl(row, function(i){
           (1+exp(beta[i]-size[i]))^(-alpha[i])
         })
  ) %>%
  group_by(size) %>%
  summarise(p2.5=quantile(prob, 0.025),
            p25=quantile(prob, 0.25),
            p75=quantile(prob, 0.75),
            p97.5=quantile(prob, 0.975),
            p50=median(prob),
            min=min(prob),
            max=max(prob))

HeV_posterior <- tibble(alpha=qunif(p=lhs[,1], min=0, max=10),
                        beta=qunif(p=lhs[,2], min=-20, max=20)) %>%
  mutate(nloglik = purrr::map2_dbl(alpha, beta, function(a, b){
    fit_observation(observed = HeV_counts$observed,
                    expected = HeV_counts$density, 
                    alpha=a, beta=b)
  }),
  lik = exp(-nloglik),
  weight = lik/max(lik),
  id=1:nrow(.)
  ) 

# resample
HeV_posterior <- HeV_posterior[sample(HeV_posterior$id, size=nsim, 
                                      prob=HeV_posterior$weight, replace=TRUE), ] %>%
  expand_grid(size=1:20) %>%
  mutate(row=1:nrow(.),
         prob = purrr::map_dbl(row, function(i){
           (1+exp(beta[i]-size[i]))^(-alpha[i])
         })
  ) %>%
  group_by(size) %>%
  summarise(p2.5=quantile(prob, 0.025),
            p25=quantile(prob, 0.25),
            p75=quantile(prob, 0.75),
            p97.5=quantile(prob, 0.975),
            p50=median(prob),
            min=min(prob),
            max=max(prob))
print(HeV_posterior[1,])

g3 <- ggplot(lik_surface_HeV)+
  geom_raster(aes(x=alpha, y=beta, fill=weight))+
  theme_classic()+
  scale_x_continuous(expression(alpha))+
  scale_y_continuous(expression(beta))+
  scale_fill_viridis_c('relative likelihood', option='A') +
  ggtitle('Hendra virus')

g4 <- ggplot(sample_dists_HeV_vis) +
  geom_segment(aes(x=size, xend=size, group=size, y=min, yend=1), 
               alpha=0.5, lwd=1) +
  geom_segment(aes(x=size, xend=size, group=size, y=p2.5, yend=1), 
               alpha=0.5, lwd=2) +
  geom_segment(aes(x=size, xend=size, group=size, y=p75, yend=p25), 
               alpha=0.5, lwd=4) +
  geom_point(aes(x=size, y=p50), size=3) +
  theme_classic()+
  scale_x_continuous('number of cases')+
  scale_y_continuous('observation probability', limits=c(0,1.1))

cowplot::plot_grid(g1, g2, g3, g4, nrow=2, labels=c('A','B','C','D'))
ggsave(file.path('figures','results.pdf'), width=8, height=5)
