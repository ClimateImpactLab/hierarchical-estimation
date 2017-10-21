Demeaning logic and time-varying covariates
-------------------------------------------

It is possible to factor out regressors that are not of interest from
only the predictor variables and not from the covariate variables.
This relies on covariates not varying between ADM2 regions within the
ADM1 regions at which fixed effects are performed, but at most varying
in time.

Suppose that we want to factor out a regressor :math:`w_{it}`, where
:math:`i` indexes regions and :math:`t` indexes time.  This is usually
done by regressing each term on it, as follows:

  :math:`x_{kit} exp(\sum_l \gamma_{kl} z_{ljt}) = \beta w_{it} +
	\epsilon_{it}`

The residuals are a linear combination of the left-hand-side
variables, according to the annihilator matrix:

  :math:`\epsilon_{it} = \sum_{hs} m_{it,hs} x_{khs} exp(\sum_l
	\gamma_{kl} z_{ljs})`

where :math:`m_{it,hs}` is an element from the annihilator matrix.
Note that :`math:`z_{ljs}` remains indexed by :math:`j` because
separate coefficients exist at ADM1 level.

If :math:`z_{ljt}` does not vary by :math:`t` and a different
coefficient exists for each region :math:`j`, this means

  :math:`\epsilon_{it} = exp(\sum_l \gamma_{kl} z_{lj}) \sum_{it} m_{it} x_{kit}`

where the term in the sum is the result of factoring out of
:math:`x_{kit}` alone.

If :math:`z_{ljt}` does vary by `t`, we can still reorganize the
expression above as

  :math:`\epsilon_{it} = \sum_t exp(\sum_l \gamma_{kl} z_{ljt}) \sum_i m_{it} x_{kit}`

The inner sum can be pre-computed, saving considerable computation
time.

