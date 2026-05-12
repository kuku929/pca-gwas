import runner
from configure import Population
from configure import config

def main():
    # generates config.yaml
    config.homogenous_population(Population.African, nsamples=100)
    config.set_output_dir(Population.African)
    config.write()
    # generates genotype information
    runner.generate_geno(nthreads=3);
    # generates plots
    runner.validate()

if __name__ == "__main__":
    main()
    # save config inside test dir for reference
    config.save_config()
