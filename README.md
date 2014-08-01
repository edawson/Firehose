Firehose_Preprocessing
======================

A collection of scripts for munging genomic and epigenetic data from
the Broad Firehose.

## Description and Motivation
The Broad Institute of MIT and Harvard provides near-monthly releases of 
pre-processed data from The Cancer Genome Atlas. Firehose remains the 
most convenient way to programatically download data from TCGA, but often 
additional preprocessing is needed to make analysis using R or python
convenient. This package provides a set of scripts for processing, cleaning,
and transforming data from the Broad Firehose site to make it easier to perform
analyses.


## Dependencies
To run any scripts in this package you'll need to install a few dependencies. First,
make sure you're running a recent version of Python (preferably from the 2.7 generation.
Next, you'll need to install a few commonly-used Python libraries for data analysis:

1. NumPy
2. SciPy
3. PANDAS
4. Matplotlib

On Ubuntu, these may in turn have some other dependencies such as:

1. python-dev
2. build-essential

The easiest way to install all of these is with a combination of APT and pip:

**NB: Always be careful and responsible using sudo**

```sudo apt-get install python-dev build-essential python-pip```

```sudo pip install numpy```

```sudo pip install scipy```

```sudo pip install pandas```

```sudo pip install matplotlib```

## Installation
Now you're ready to install the Firehose package:

```git clone https://github.com/edawson/Firehose.git```

Right away you can run scripts out of the Firehose directory, but if you'd like
to link them into /usr/bin/ and run them anywhere you can do this:

```sudo ln -s "<PATH/TO/DESIRED/SCRIPT/script_name.py" "/usr/bin/script_name```


## Usage
Usage should be consistent across the entire package (with one caveat we'll get to later):

```python name_of_script.py -i <FILE_FOR_ANALYSIS>```

That's it! Currently there are no options for fine-tuning processing because these scripts try
to stay as minimalistic as possible. As time goes on usage may be subject to change, but
the above usage will **always** perform the standard transformations.

## Contacting the Author
If you use these scripts please make sure to credit where they came from:

>Eric T. Dawson</br>
>github.com/edawson/Firehose<br/>
>erictdawson.com

If you would like to contact the author for comments, criticism, or praise, he can be reached
through his website above.

## Acknowledgements
These scripts were developed during the author's stint as a student intern at The Ontario
Institute for Cancer Research. They are a (nearly) direct port of scripts originally produced
by Dr. Guanming Wu.
