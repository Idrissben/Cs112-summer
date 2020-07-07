from maccabee.constants import Constants
from maccabee.data_sources.data_source_builders import build_random_normal_datasource
from maccabee.data_generation.data_generating_process_sampler import DataGeneratingProcessSampler
from maccabee.data_generation.data_generating_process import SampledDataGeneratingProcess

from maccabee.parameters import build_parameters_from_axis_levels
from maccabee.parameters import build_default_parameters

# pprint is used to examine the classes used throughout - uncomment relevant lines to see
# from pprint import pprint

HIGH, MEDIUM, LOW = Constants.AxisLevels.HIGH, Constants.AxisLevels.MEDIUM, Constants.AxisLevels.LOW

# iterate over all possible settings of a parameter
for i in [HIGH, MEDIUM, LOW]:
    dgp_params = build_parameters_from_axis_levels({
        Constants.AxisNames.OVERLAP: i,                   # choose desirable parameter to manipulate here
        Constants.AxisNames.PERCENT_TREATED: LOW          # LOW defaults to 35% - use the next line to specify
    })
    
    # dgp_params.set_parameter("TARGET_PROPENSITY_SCORE", 0.8)
    
    # pprint(vars(dgp_params))
    
    normal_data_source = build_random_normal_datasource(
        n_covars=5,
        n_observations=5000)
    
    # pprint(vars(normal_data_source))


    dgp_sampler = DataGeneratingProcessSampler(
        dgp_class=SampledDataGeneratingProcess,
        parameters=dgp_params,
        data_source=normal_data_source,
        dgp_kwargs={})
    
    sampled_dgp = dgp_sampler.sample_dgp()
    # pprint(vars(sampled_dgp))
    sampled_dataset = sampled_dgp.generate_dataset()
    
    # use this line to display the generated dataset in python
    # sampled_dataset.observed_data 
    # use this line to save the generated dataset to file; use string formatting if you're saving multiple files in a single loop
    # sampled_dataset.observed_data.to_csv("data_overlap_{}.csv".format(i))
