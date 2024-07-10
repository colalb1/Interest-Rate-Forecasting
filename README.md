# Interest Rate Forecasting

Implementation of common and state-of-the-art interest rate models in C++. This project aims to act as an implementation and mathematical reference for interest rate models that build upon the basics (Vasicek, CIR, etc.) via more complex dynamics.

## Files

## Models

If you are unfamiliar with the classic interest rate models, I recommend reading [this](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) book before continuing. This text was used for mathematical reference of the classic models.

### Classic Models
**Vasicek:**

The Vasicek model states that the instantaneous interest rate is derived from the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dr(t) = \kappa(\theta - r(t))dt + \sigma dW(t)$$

where $r(0) = r_0$. This is an Ornstein-Uhlenbeck process with constant coefficients under a risk-neutral measure, meaning this describes a [Gauss-Markov](https://en.wikipedia.org/wiki/Gauss%E2%80%93Markov_process) process (stationary process based on the normal distribution) of a variable that tends to revert to some mean ($\theta$) over time given assets follow [risk-neutral pricing](https://en.wikipedia.org/wiki/Risk-neutral_measure) principles. The variables are assumed to be constant; this is unrealistic, but this is a basic model. The variables are defined as follows:

* $\kappa > 0$ is the reversion rate (how fast the process reverts to the mean)
* $\theta$ is the long-term mean to which the process reverts
* $\sigma > 0$ is the volatility coefficient
* $W(t)$ is a standard Brownian motion
* $r(t)$ is the short-term interest rate that one aims to solve for

The $\sigma dW(t)$ term represents the random fluctuations about the walk, and the $\kappa(\theta - r(t))$ term represents the drift term that reverts the interest rate to the given mean ($\theta$).

One may solve the differential equation conditional on [filtration](https://en.wikipedia.org/wiki/Filtration_(probability_theory)) $\mathcal{F}(s)$ where $s\leq t$ and obtain the following analytic terms for the expectation and variance of the walk:

$$\mathbb{E}(r(t) | \mathcal{F}(s)) = r(s)\exp(-\kappa (t - s)) + \theta (1 - \exp(-\kappa (t - s)))$$

$$Var(r(t) | \mathcal{F}(s)) = \frac{\sigma ^ 2}{2\kappa}\left(1 - \exp(-2\kappa (t - s)\right)$$

Given the simplicity of the model, one may derive the price of a pure-discount bond given current time $t$ and expiration date $T$:

$$P(t, T) = \exp\left[A(t, T) - B(t, T)r(t)\right]$$

where

$$B(t, T) = \frac{1}{k}\left(1 - \exp(1 - \kappa (T - t))\right)$$

$$A(t, T) = \left(\theta - \frac{\sigma ^ 2}{2\kappa ^ 2}\right)(B(t, T) - T + t) - \frac{\sigma ^ 2}{4k}B(t, T)^2$$

I will defer further mathematical explanation to Section **3.2.1** of [the reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) for brevity.

A few main issues of the Vasicek model are that it allows for negative interest rates due to the unconstrained movement of the $\sigma dW(t)$ term, unrealistic constant volatility/reversion speed, and the assumption that the process is based on the normal distribution. The following model addresses the negative interest rates.

**Cox-Ingress-Ross (CIR):**

The CIR model addresses the issue of negative interest rates in the Vasicek model by adding a $\sqrt{r(t)}$ term to the random fluctuation term. Thus, the process is modeled by the following stochastic differential equation:

$$dr(t) = \kappa(\theta - r(t))dt + \sqrt{r(t)} \sigma dW(t)$$

This prevents negative interest rates since the steps become infinitesimally small when $r(t)\to 0$ and must revert from $0$ for the process to continue. Essentially, the steps get smaller by a factor of $\sqrt{r(t)}$ when $r(t)\to 0$ so it is impossible to reach a negative rate.

All other variables in the CIR model are defined similarly to the Vasicek model. Although slightly more complex, the CIR model is still simple enough to derive the expected rate, variance, and price of a pure-discount bond analytically:

$$\mathbb{E}(r(t) | \mathcal{F}(s)) = r(s)\exp(-\kappa (t - s)) + \theta (1 - \exp(-\kappa (t - s)))$$

$$Var(r(t) | \mathcal{F}(s)) = r(s)\frac{\sigma ^ 2}{\kappa}\exp(-\kappa (t - s))\left(1 - \exp(-\kappa (t - s)\right) + \theta\frac{\sigma ^ 2}{2\kappa}\left(1 - \exp(-\kappa (t - s))\right)^2$$

$$P(t, T) = A(t, T)\exp(-B(t, T)r(t))$$

where

$$A(t, T) = \left[\frac{2h\exp\left(\frac{(k + h)(T - t)}{2}\right)}{2h + (k + h)(\exp((T - t)h) - 1)}\right]^{\frac{2\kappa\theta}{\sigma ^ 2}}$$

$$B(t, T) = \frac{2(\exp((T - t)h) - 1)}{2h + (k + h)(\exp((T - t)h) - 1)}$$

$$h = \sqrt{\kappa ^ 2 + 2\sigma ^ 2}$$

I will defer to Section **3.2.3** of the [reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) for further explanation.


The CIR model suffers from similar issues to the Vasicek model apart from negative rates. Extreme movements influence real-life financial data to observe distributions with fatter tails than the normal distribution, implying models such as the CIR model underestimate risk. The **Stable CIR** model will address this issue by basing the movement term on a [Levy alpha-stable distribution](https://en.wikipedia.org/wiki/Stable_distribution) instead of the normal distribution.

### Modern Models
**Stable CIR:**

**$\alpha$-CIR**



Using [this](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) as reference.

## TODO:

- WRITE README FOR COMPLETED PARTS
- Improve mathematical formulation of current models
- Implement improvement of current model
- Test improved current model
- Write technical README
- Write intro and conclusion README
- Share

USE INTEL VTUNE TO FIND OPTIMAL LOOP-BLOCKING SIZE
