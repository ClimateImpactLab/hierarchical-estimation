'''
Estimates global mortality response by year*adm2*age

:author: My Name
:contact: `myname@gmail.com <mailto:myname@gmail.com>`_
:modified: 2016/08/19
:team: Mortality
:lead: Amir Jina


Input variables
~~~~~~~~~~~~~~~

global_mortality_admin2_year_3agegrp.dta
    :location: ~/Dropbox/ChicagoRAs_Reanalysis/interpolation/data/consolidated
    :description: Global baseline mortality data
    :version: GCP.Mortality.GlobalBaselineMortality.2016-08-19


Output variables
~~~~~~~~~~~~~~~~

global_6countries_response_adm0prcp_global_year_FE.dta
    :location: ~/Dropbox/GCP/MORTALITY/tables/Table_2/Panel_A/ster
    :description: 6-country response surface with year FE
    :version: ``GCP.Mortality.Global_6countries_response_adm0prcp_yearFE.2016-08-19``

    The year FE dataset uses a linear regression with global year fixed effects:

        :math:`M_{i, t}=\sum_{i\in K}{\hat{\beta}_{i, t}^{k}T_{i, t}^{k}X_{i, t}}+\delta_{t}+\varepsilon_{i, t}`

'''

def multiply_numbers(*numbers):
    '''
    Takes arguments and returns the product of these arguments

    Parameters
    ----------
    numbers : float
        Floats to multiply

    Returns
    -------
    float
        Product of all arguments

    .. math::
        \Pi_{i \in #(numbers)}{numbers_i}

    '''

    return reduce(lambda x, y: x*y, *floats)
