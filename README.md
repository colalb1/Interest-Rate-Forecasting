# Interest Rate Forecasting

Implementation of common and state-of-the-art interest rate models in C++. This project aims to act as an implementation and mathematical reference for interest rate models that build upon the basics (Vasicek, CIR, etc.) via more complex dynamics.

## Files

## Models

If you are unfamiliar with the classic interest rate models, I recommend reading [this](https://www.amazon.com/Interest-Rate-Models-Practice-Inflation/dp/3540221492) book before continuing. This text is used for reference to the classic models.

### Classic Models
**Vasicek:**

The Vasicek model states that the instantaneous interest rate is derived from the following [stochastic differential equation](https://en.wikipedia.org/wiki/Stochastic_differential_equation):

$$dr(t) = \kappa(\theta - r(t))dt + \sigma dW(t)$$

where $r(0) = r_0$.

**Cox-Ingress-Ross (CIR):**

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

USE INTEL VTUNE TO FIND OPTIMAL LOOP BLOCKING SIZE
