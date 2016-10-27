%{
Generate correlated errors within age-groups

:author: My Name
:contact: `myname@gmail.com <mailto:myname@gmail.com>`_
:modified: 2016/08/19
:team: Mortality
:lead: Amir Jina


.. math::

    \mathbb{vcv} = \left[\begin{array}[b]{ccc}C\sigma_1^2\mathbb{1}+(1-c)\sigma_1^2\mathbb{e} & 0 & 0 \\0 & C\sigma_1^2\mathbb{1}+(1-c)\sigma_1^2\mathbb{e} & 0 \\0 & 0 & C\sigma_1^2\mathbb{1}+(1-c)\sigma_1^2\mathbb{e}\end{array}\right]

%}

clear
clc

%st dev of error by age-group
sd31 = 250;
sd32 = 250;
sd33 = 1000;
%within age group correlation
co = 0.9;

%define vcv matrix
cov31 = co*sd31^2*ones(4590,4590)+(1-co)*sd31^2*eye(4590);
cov32 = co*sd32^2*ones(4590,4590)+(1-co)*sd32^2*eye(4590);
cov33 = co*sd33^2*ones(4590,4590)+(1-co)*sd33^2*eye(4590);
vcv = blkdiag(cov31,cov32,cov33);

%generate 5000 draws of 13770 errors
mu = zeros(1,13770);
rng default  % For reproducibility
draws = mvnrnd(mu,vcv,5000);
errors = transpose(draws);

corrcoef(draws(:,1:2))
corrcoef(draws(:,4590:4592))

