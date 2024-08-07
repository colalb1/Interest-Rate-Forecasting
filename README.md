# Interest Rate Forecasting

Implementation of common interest rate models and (a few of) their extensions in C++. This project aims to act as an implementation and mathematical reference for interest rate models that build upon the basics (Vasicek, CIR, etc.) via more complex dynamics.

## Files

The [src](https://github.com/colalb1/Interest-Rate-Forecasting/tree/main/src) folder contains the cpp files with the models.

[*Vasicek-CIR.cpp*](https://github.com/colalb1/Interest-Rate-Forecasting/blob/main/src/Vasicek-CIR.cpp) contains the basic Vasicek, exponential Vasicek, and CIR implementations.

[*General-CIR.hpp*](https://github.com/colalb1/Interest-Rate-Forecasting/blob/main/src/General-CIR.cpp) contains the Stable CIR and $\alpha$-CIR implementations.

[*BK.cpp*](https://github.com/colalb1/Interest-Rate-Forecasting/blob/main/src/BK.cpp) contains the Black-Karasinski implementation.

[*HW.cpp*](https://github.com/colalb1/Interest-Rate-Forecasting/blob/main/src/HW.cpp) contains the One-Factor and Two-Factor Hull-White models.

## Models

If you are unfamiliar with the classic interest rate models, I recommend reading [this](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) book; this text was used for mathematical reference of the classic models. A basic description of the classic models is also given below. If you are familiar with them, skip to the **Modern Models** subsection.

The [Euler-Maruyama method](https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method) was used for all numerical integration.

### Classic Models
**Vasicek:**

The Vasicek model states that the instantaneous interest rate is derived from the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dr(t) = \kappa(\theta - r(t))dt + \sigma dW(t)$$

where $r(0) = r_0$. This is an Ornstein-Uhlenbeck process with constant coefficients under a risk-neutral measure, meaning this describes a [Gauss-Markov](https://en.wikipedia.org/wiki/Gauss%E2%80%93Markov_process) process (stationary process based on the normal distribution) of a variable that tends to revert to some mean $\theta$ over time given assets follow [risk-neutral pricing](https://en.wikipedia.org/wiki/Risk-neutral_measure) principles. The variables are assumed to be constant; this is unrealistic, but this is a basic model. The variables are defined as follows:

* $\kappa > 0$ is the reversion rate (how fast the process reverts to the mean)
* $\theta$ is the long-term mean to which the process reverts
* $\sigma > 0$ is the volatility coefficient
* $W(t)$ is a standard Brownian motion
* $r(t)$ is the short-term interest rate that one aims to solve for

The $\sigma dW(t)$ term represents the random fluctuations about the walk, and the $\kappa(\theta - r(t))$ term represents the drift term that reverts the interest rate to the given mean ($\theta$).

One may solve the differential equation conditional on [filtration](https://en.wikipedia.org/wiki/Filtration_(probability_theory)) $\mathcal{F}(s)$ where $s\leq t$ and obtain the following analytic terms for the expectation and variance of the walk:

$$\mathbb{E}(r(t) | \mathcal{F}(s)) = r(s)e^{-\kappa(t - s)} + \theta(1 - e^{-\kappa(t - s)})$$

$$Var(r(t) | \mathcal{F}(s)) = \frac{\sigma ^ 2}{2\kappa}\left(1 - e^{-2\kappa(t - s)}\right)$$

Given the simplicity of the model, one may derive the price of a pure-discount bond given current time $t$, expiration date $T$, and time to expiry $\tau = T - t$:

$$P(t, T) = \exp\left[A(t, T) - B(t, T)r(t)\right]$$

where

$$B(t, T) = \frac{1}{k}\left(1 - e^{-\kappa\tau}\right)$$

$$A(t, T) = \left(\theta - \frac{\sigma ^ 2}{2\kappa ^ 2}\right)(B(t, T) - \tau) - \frac{\sigma ^ 2}{4k}B(t, T)^2$$

I will defer further mathematical explanation to Section 3.2.1 of [the reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) for brevity.

A few issues of the Vasicek model are that it allows for negative interest rates due to the unconstrained movement of the $\sigma dW(t)$ term, unrealistic constant volatility/reversion speed, and the assumption that the process is based on the normal distribution. The CIR model addresses the negative interest rate issue.

**Cox-Ingress-Ross (CIR):**

The CIR model addresses the issue of negative interest rates in the Vasicek model by multiplying the random fluctuation term by $\sqrt{r(t)}$. Thus, the process is modeled by the following stochastic differential equation:

$$dr(t) = \kappa(\theta - r(t))dt + \sigma \sqrt{r(t)} dW(t)$$

This prevents negative interest rates since the steps become infinitesimally small when $r(t)\to 0$ and must revert from $0$ for the process to continue.

All other variables in the CIR model are defined similarly to the Vasicek model. Although slightly more complex, the CIR model is still simple enough to derive the expected rate, variance, and price of a pure-discount bond analytically:

$$\mathbb{E}(r(t) | \mathcal{F}(s)) = r(s)e^{-\kappa (t - s)} + \theta (1 - e^{-\kappa (t - s)})$$

$$Var(r(t) | \mathcal{F}(s)) = r(s)\frac{\sigma ^ 2}{\kappa}e^{-\kappa (t - s)}\left(1 - e^{-\kappa (t - s)}\right) + \theta\frac{\sigma ^ 2}{2\kappa}\left(1 - e^{-\kappa (t - s)}\right)^2$$

$$P(t, T) = A(t, T)\exp(-B(t, T)r(t))$$

where

$$A(t, T) = \left[\frac{2he^{\frac{(k + h)\tau}{2}}}{2h + (k + h)(e^{h\tau} - 1)}\right]^{\frac{2\kappa\theta}{\sigma ^ 2}}$$

$$B(t, T) = \frac{2(e^{h\tau} - 1)}{2h + (k + h)(e^{h\tau} - 1)}$$

$$h = \sqrt{\kappa ^ 2 + 2\sigma ^ 2}$$

I will defer to Section 3.2.3 of the [reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) for further explanation.


The CIR model suffers from issues similar to those of the Vasicek model, apart from negative rates. Extreme movements influence real-life financial data to observe distributions with fatter tails than the normal distribution, implying models such as the CIR model underestimate risk. The Stable CIR model will address this issue by basing the movement term on a [Levy alpha-stable distribution](https://en.wikipedia.org/wiki/Stable_distribution) instead of the normal distribution.

**Exponential Vasicek (EV):**

The EV model assumes the short rate follows a [lognormal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) instead of a normal distribution as the Vasicek model does. Suppose that $\log(r(t))$ follows an Ornstein-Uhlenbeck process $y(t)$. By definition, the Vasicek SDE describes $y(t)$ as follows:

$$dy(t) = (\theta - \kappa y(t))dt + \sigma dW(t)$$

By [Ito's Lemma](https://en.wikipedia.org/wiki/It%C3%B4%27s_lemma), the SDE for $r(t)$ is the following:

$$dr(t) = r(t)\left(\theta + \frac{\sigma ^ 2}{2} - \kappa\log(r(t))\right)dt + \sigma r(t)dW(t)$$

The explicit definition of the expectation and variance of $r(t)$ under this model will be left to page 71 of the [reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492); the analog under filtration $\mathcal{F}(s)$ as $t\to\infty$ will be given instead. This essentially gives the expected value and variance for long-term predictions.

$$\lim_{t\to\infty}\mathbb{E}(r(t) | \mathcal{F}(s)) = \exp\left(\frac{\theta}{\kappa} + \frac{\sigma ^ 2}{4\kappa}\right)$$

$$\lim_{t\to\infty}\mathbb{Var}(r(t) | \mathcal{F}(s)) = \exp\left(\frac{2\theta}{\kappa} + \frac{\sigma ^ 2}{2\kappa}\right)\left[\exp\left(\frac{\sigma ^ 2}{2\kappa} - 1\right)\right]$$

There are no explicit formulas for pure-discount bonds; they can be solved numerically instead.

The EV model suffers from an infinite expected value of the money-market account regardless of maturity as $r(t)$ is assumed to be lognormal. This is a realism issue for all extensions based on the EV model.


**Black-Karasinski (BK)**

The BK model is an extension of the EV model where the mean $\theta$ is a function that depends on time. Thus, the short rate is assumed to follow a lognormal distribution, and its SDE is given by the following:

$$dr(t) = r(t)\left(\theta(t) + \frac{\sigma ^ 2}{2} - \kappa\log(r(t))\right)dt + \sigma r(t)dW(t)$$

No analytic formulas for discount bonds or option bonds exist, similar to the EV model, although they may be simulated with a [trinomial tree](https://en.wikipedia.org/wiki/Trinomial_tree).

Similar to the EV model, the expected money-market account value is infinite regardless of maturity since it is assumed that $r(t)$ follows a lognormal distribution.

**One-Factor Hull-White (HW1)**

The HW1 model is an extension of the Vasicek model as it allows the long-term mean $\theta$ and volatility $\sigma$ to vary with time. This allows for more precise calibration of the current term structure and can accommodate various [yield curve](https://www.investopedia.com/terms/y/yieldcurve.asp) shapes, adding a layer of realism.

This model assumes $r(t)$ follows a normal distribution, and its SDE follows as such:

$$dr(t) = (\theta(t) - \kappa r(t))dt + \sigma(t)dW(t)$$

Although more complex than the Vasicek model, HW1 is still too simple to capture the intricacies of the entire yield curve given interest rate movements typically depend on multiple sources of risk (multiple Brownian motions). Extensions of the HW1 model, such as the Two-Factor Hull-White model or Multi-Factor Hull-White models were developed to address this issue. The following section will detail the Two-Factor Hull-White model. 

**Two-Factor Hull-White (HW2)**

As mentioned in the previous section, HW2 extends HW1 by adding a second risk (stochastic factor, Brownian motion in this case) to provide additional flexibility to the term structure. The SDEs that define this model are the following:

$$dr(t) = (\theta(t) - \kappa_1r(t) + x(t))dt + \sigma_1(t)dW_1(t)$$
$$dx(t) = -\kappa_2x(t)dt + \sigma_2(t)dW_2(t)$$
$$dW_1dW_2 = \rho dt$$

$r(t)$ is the short rate and $x(t)$ is an additional stochastic process that influences the $r(t)$ dynamics. Similar to the unindexed definitions, $\kappa_1, \kappa_2$ represent the mean reversion rates of their respective equations while $\sigma_1, \sigma_2$ are defined similarly for volatility. $\rho$ is the correlation coefficient between the Brownian motions. The correlation term allows for a more thorough assessment of joint risk and improves hedging strategies via completeness of information.

HW2 is equivalent to the [G2++ model](https://www.cequra.uni-muenchen.de/download/cequra_wp21.pdf); G2++ is typically used in practice as its implementation is simpler than HW2. I will leave further mathematical explanation to section 4.2 of the [reference book](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492).

### Modern Models

The mathematics of the Stable CIR and $\alpha$-CIR models are derived primarily from [this](https://arxiv.org/abs/2402.07503) paper. 

**Stable CIR:**

As mentioned in the CIR subsection, the Stable CIR model aims to enhance the classic CIR model by basing the motion on fatter-tailed distributions as real data tends to show larger aberrations than that of idealized models using the normal distribution via Brownian motion. [This paper](https://arxiv.org/abs/1301.3243) gives more mathematical background to WHY we want a fatter-tailed distribution; I encourage you to read this for more background, but I will be omitting the explanation and will move on to the HOW of the problem.

First, I must define the Levy alpha-stable distribution class. Four parameters define said class of distributions: the stability parameter $\alpha\in(0, 2]$, skewness parameter $\beta\in[-1, 1]$, scale parameter $\sigma\in(0, \infty)$, and the location parameter $\mu\in(-\infty, \infty)$. Put simply, $\mu$ is the mean/center of the distribution, $\sigma$ is the statistical dispersion factor (variance is a common example), $\beta$ measures asymmetry ($\beta = 0$ implies symmetry while $\beta = -1, 1$ imply left and right skewness, respectively), and $\alpha$ measures tail heaviness (lower $\alpha\implies$ heavier tails). $\alpha$ is the parameter of main interest. The following facts about the $\alpha$-stable distribution are relevant to this context.

* Distributions with $\alpha = 2$ are normal
* Distributions with $\alpha = 1$ are [Cauchy](https://en.wikipedia.org/wiki/Cauchy_distribution)
* Distributions with $\alpha < 2$ have undefined second and higher moments (variance and higher moments)
* Distributions with $\alpha \leq 1$ have undefined first and higher moments (mean and higher moments)

Given this information, one must limit the distributions to $\alpha\in(1, 2]$ so the drift term(s) is/are well-defined (the drift term needs a mean). There is no general analytic form of the probability distribution function, but its [characteristic function](https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)) is defined as follows (taken from the [alpha-stable distribution Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)):

$$\phi(t, \alpha, \beta, \sigma, \mu) = \exp\left(it\mu - |\sigma t|^{\alpha}(1 - i\beta sgn(t)\Phi)\right)$$

where 

$$\Phi = \begin{cases}
\tan\left(\frac{\pi\alpha}{2}\right) \text{, } \alpha \neq 1\\
-\frac{2}{\pi}\log|t| \text{, } \alpha = 1
\end{cases}$$

The characteristic function is not particularly useful for this application; it is provided for context and your general enlightenment.

Thus, the [inverse transform method](https://en.wikipedia.org/wiki/Inverse_transform_sampling) cannot be used since the cumulative distribution function (cdf) cannot be derived. Instead, the inverse cdf can be estimated with the stochastic model from [this paper](https://www.sciencedirect.com/science/article/pii/0167715295001131). Said paper did not hypothesize the model but proved its correctness. The model is defined as follows (VERY nested, but it is relatively straightforward):

$$Y = \begin{cases}
\mu + \sigma X \text{, } \alpha\neq 1,\\
\mu + \sigma X + \frac{2}{\pi}\beta\sigma\log(\sigma) \text{, }\alpha = 1
\end{cases}$$

given 

$$X = \begin{cases}
S_{\alpha, \beta} \left(\frac{\sin(\alpha(V + B_{\alpha, \beta}))}{(\cos(V))^{1 / \alpha}}\right) \left(\frac{\cos(V - \alpha(V + B_{\alpha, \beta}))}{W}\right)^{(1 - \alpha) / \alpha} \text{,  }\alpha\neq 1\\
\frac{2}{\pi}\left[\left(\frac{\pi}{2} + \beta V\right)\tan(V) - \beta\log\left(\frac{W\cos(V)}{\pi / 2 + \beta V}\right)\right] \text{,  }\alpha = 1
\end{cases}$$

where 

$$B_{\alpha, \beta} = \frac{\arctan(\beta\tan(\frac{\pi\alpha}{2}))}{\alpha}$$

and 

$$S_{\alpha, \beta} = \left(1 + \beta ^ 2\tan^2\left(\frac{\pi\alpha}{2}\right)\right)^{1 / (2\alpha)}$$

given random variables

$$\begin{cases}
V\sim Unif\left(-\frac{\pi}{2}, \frac{\pi}{2}\right)\\
W\sim Exp(1)\end{cases}$$

The uniform and exponential distribution generators in C++ simulated these random variables. This is known as the **Chamber-Mallows-Stuck** method.

Most of the legwork is done defining the distribution; the model is straightforward in comparison. The stochastic differential equation is defined as follows:

$$dr(t) = (\kappa r(t-) + \theta)dt + (\sigma r(t-))^{1 / \alpha} dZ^{\alpha}(t)$$

where $Z^{\alpha}(t)$ is an alpha-stable process and $r(t-)$ is the rate at the previous timestep.

This is essentially the same as the CIR model except the random walk is based on an alpha-stable distribution instead of a normal distribution. The $\alpha$-CIR model will generalize this by adding additional $(\sigma r(t-))^{1 / \alpha} dZ^{\alpha}(t)$ terms to more accurately reflect real-life market conditions via multiple independent rate-dependent risks.

**$\alpha$-CIR**

The SDE of this process is defined as follows:

$$dr(t) = (\kappa r(t-) + \theta)dt + \sum_{i = 1}^g (\sigma_ir(t-))^{1 / \alpha_i}dZ_i^{\alpha_i}(t)$$

given indices $1 < \alpha_1 < \dots < \alpha_g \leq 2$, independent stable processes $Z_i^{\alpha_i}$, and rate-dependent volatilities $\sigma_i$ for $i = 1, \dots g$.

The summation term represents multiple independent sources of stochastic noise that scale based on the current interest rate, each process with a unique distribution. The $(\sigma_ir(t-))^{1 / \alpha_i}$ term is state-dependent volatility that provides additional realism via [volatility clustering](https://en.wikipedia.org/wiki/Volatility_clustering#:~:text=In%20finance%2C%20volatility%20clustering%20refers,be%20followed%20by%20small%20changes.%22). Each $(\sigma_ir(t-))^{1 / \alpha_i}dZ_i^{\alpha_i}(t)$ term represents a different risk, allowing for varying degrees of tail behavior and jump densities for a more nuanced representation of market risks.

In short, one may control the $\alpha_i$ and $\sigma_i$ terms for each risk to represent market conditions accurately at a given interest rate. If you are interested in more of the statistical-theoretic nuance (classification of generating equations, canonical representation, moments of the rates, etc.) of this model, read Section 3 of [the paper](https://arxiv.org/abs/2402.07503).

This model may be improved by removing the assumption of independence between risks as real-life risks are seldom independent.

## Programmatic Optimizations

**Loop-blocking**

Also known as loop-tiling, this improves computation speed by breaking down large blocks/loops into smaller blocks/tiles that fit into various levels of cache more efficiently. Smaller blocks of data allow the data to stay in the smaller, faster caches (L1 or L2) and reduce cache misses since there are fewer evictions. This allows for multi-level cache utilization that can be tuned using programs such as [Intel VTune](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html). One may adjust the block size constant for optimal efficiency.


## Takeaways

The goal of this project was to refresh and expand my knowledge of interest rate models, create reference material for those interested in the basics of modeling interest rate SDEs, and further improve my C++ knowledge/OOP organization.

Although this repo mainly contains well-studied interest rate SDEs, this will act as the basis for more nuanced implementations moving forward. I will come back to this project later to add additional modern academic techniques; admittedly, I did not get to as much of this content as I would have liked. Plans for the future of this (or a tangential) project include focusing on extensions of multi-factor models such as the [Longstaff-Schwartz](https://www.soa.org/493589/globalassets/assets/library/newsletters/risks-and-rewards/2002/october/rar-2002-iss40-smith.pdf), [Chen](https://en.wikipedia.org/wiki/Chen_model), [Heath-Jarrow-Merton family](https://en.wikipedia.org/wiki/Heath%E2%80%93Jarrow%E2%80%93Morton_framework) (HJM), and [LIBOR](https://en.wikipedia.org/wiki/LIBOR_market_model) models (a commonly used subset of the HJM framework).

A significant amount of time was spent reading C++ documentation to ensure the code reached an adequate standard of efficiency and formatting. The standards/guidelines I followed were the [CppCoreGuidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) and [muOS++](https://micro-os-plus.github.io/develop/sutter-101/). [Effective C++](https://www.amazon.com/Effective-Specific-Improve-Programs-Designs/dp/0321334876) was used supplementally for code design. That being said, there is room for improvement. If you are experienced with C++ and would like to leave some comments or do a code review, please contact me on [LinkedIn](https://www.linkedin.com/in/colin-alberts/).


## Side Note
The descriptions came out much longer than I intended, but I figured I would make it detailed since my only other life responsibility while completing the bulk of this project was playing Elden Ring. At the time of this writing, I have a level 200+ Dex build using the [Backhand Blade](https://eldenring.wiki.fextralife.com/Backhand+Blade). Edit: I had to change back to my heavy affinity [Greatsword](https://eldenring.wiki.fextralife.com/Greatsword) build for Consort Radahn; it was quite difficult.

If you are interested in more of the project details or my Elden Ring build, please contact me via [LinkedIn](https://www.linkedin.com/in/colin-alberts/).
