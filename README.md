# Random_samples_bootstrapping_clustering

## **What is this repository for?**

In this script, I create 100 random samples following the bootstrapping approach described in Ling et al. 1986. I create pseudo-data samples by randomly drawing data from our parent sample, allowing for repetitions. This is very useful when further data samples are needed to perform specific tests or when poissonian error bars may understimate the uncertainties of the data. We also calculate the clustering of a data sample with the K-estimator of Adelberger et al. 2005, comparing data pairs along two different radial separation bins at fixed transverse separations. I finally compare the bootstrapping and poissonian error bars.

## **Installing Random_samples_bootstrapping_clustering**

No installation is needed. Simply cloning the GitHub repository and importing the script is enough. Hence, 

```
    git clone https://github.com/YohanaHerrero/Random_samples_bootstrapping_clustering.git
```

The code is written in Python and uses a range of default packages included in standard installations of Python:

### **Standard packages:**

- numpy  
- matplotlib
- math
- glob
- itertools
- sklearn

### **Special packages:**

- astropy 

After adding the Random_samples_bootstrapping_clustering directory to the PYTHONPATH or changing location to the Random_samples_bootstrapping_clustering directory, the repository can be imported in python with:

```
    import Random_samples_bootstrapping_clustering
```

Besides the python packages, you will also need a fits table containing the following columns:

- RA in degrees
- DEC in degrees
- Z 

Instead, the three parameters above can be replaced by various kind of data, according to specific needs.
