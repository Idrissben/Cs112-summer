from maccabee.constants import Constants
from maccabee.data_sources.data_source_builders import build_random_normal_datasource
from maccabee.data_generation.data_generating_process_sampler import DataGeneratingProcessSampler
from maccabee.data_generation.data_generating_process import SampledDataGeneratingProcess

from maccabee.parameters import build_parameters_from_axis_levels
from maccabee.parameters import build_default_parameters

# pprint is used to examine the classes used throughout - uncomment relevant lines to see
# from pprint import pprint

HIGH, MEDIUM, LOW = Constants.AxisLevels.HIGH, Constants.AxisLevels.MEDIUM, Constants.AxisLevels.LOW

for i in [HIGH, MEDIUM, LOW]:
    dgp_params = build_parameters_from_axis_levels({
        Constants.AxisNames.TE_HETEROGENEITY: i
    })
    
    dgp_params.set_parameter("TARGET_PROPENSITY_SCORE", 0.05)
    
    # pprint(vars(dgp_params))
    
    normal_data_source = build_random_normal_datasource(
        n_covars=3,
        n_observations=1000)
    # pprint(vars(normal_data_source))


    dgp_sampler = DataGeneratingProcessSampler(
        dgp_class=SampledDataGeneratingProcess,
        parameters=dgp_params,
        data_source=normal_data_source,
        dgp_kwargs={})
    
    sampled_dgp = dgp_sampler.sample_dgp()
    # pprint(vars(sampled_dgp))
    sampled_dataset = sampled_dgp.generate_dataset()
    
    sampled_dataset.observed_data.to_csv("trt_hetero_{}.csv".format(i))
