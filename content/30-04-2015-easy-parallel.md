Title: Easy Parallel Computing with IPython
Date: 30-04-2015
Authors: Andy Mellor
Slug: easy-parallel
Tags: IPython, Parallel, Fractals
Category: Python
Status: Published

<!-- PELICAN_BEGIN_SUMMARY -->

After speaking to a number of academics involved with Monte Carlo simulations, I realised many of them resort to creating bash scripts to 
run a number of python instances of simulation before writing another script to piece the data back together. It turns out this is incredibly simple 
to do in IPython, and in particular in an IPython Notebook. We take a look at how to run computations in parallel as well as giving a use-case in the creation of Julia fractals.

<!-- PELICAN_END_SUMMARY -->

Most modern CPUs have 2 or more cores, allowing for everyone to take full advantage of parallel computing. The IPython implementation here is probably most useful for when
you have some computationally intensive program which is not called on a regular basis. Writing a full parallel implementation can be time consuming (depending on the problem and the 
experience of the programmer), however as we'll see, writing an implementation in IPython can be quick and painless. Of course more tailored parallel implementations require more effort 
should use a more established MPI such as OpenMPI. Thankfully this can be done within the IPython framework.

To begin, we will set up the engines that we will use for our computation. 

{% notebook parallel.ipynb %}


