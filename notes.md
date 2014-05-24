Extra notes on `ridgeds`
========================


*What am I trying to show here?* Three things:

 1. what is the effect of covariates that are highly correlated with distance on estimate of detectability?
 2. What about data that isn't correlated at all?
 3. what about covariates that are highly correlated with each other?

Last one is probably most interesting to most people.

Thoughts: 
  * How do we simulate these situations?
  * Not really a problem if we don't select the variable by AIC?




## 1. Covariates correlated with distance

What would be happening mathematically? Covarite is a mix of the true covariate and the observed?

What is a realistic model of this issue?

Simulate "visability" or "sighting conditions"-type data.

Some "transects" (doesn't matter, since the detection function doesn't take any notice) with "good" sightability, some with variable sightability. Multinomial w. higher weight on zero?


Becker and Quang (2009) JABES (p. 211) has a description of search distance -- seems like a truncation, maybe generate observations using different truncations?

## 2. junk data

May as well do this as uniform(0,1) as we can always normalise the covariate values? Maybe uniform(-1,1) ?

## 3. correlated covariates

Which model comes out AIC best?

Compare $P_a$:
  * bias?
  * variance?

Compare parameter estimates

