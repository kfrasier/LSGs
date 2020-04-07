Step 1: run LTSAdailySpectra_step1.m

line 24: modify the path to point to HARP data summary CSV.
line 27: modfy to point to your TFs folder.

The code will TRY to use this info to suggest the correct transfer function file, but pay attention!! 
It's not perfect, and if you are testing new TFs, it might point to a spot you didn't intend.

Run this on all LTSAs from a given deployment (select multiple files when prompted)
This is the slowest step.


Step 2: run plotDailyAveSpectra_step2.m

I made a _df20 version of this because I needed to modify number of frequency bins, but it would be possible to just make the original adaptable.

This step requires a _pltparams.txt file. See example in the folder.
Variables include:  
     ifile = modify to correct input file name (produced by step 1)
     ipath = modify to correct input path name
     navepd = number of averages per day. Depends on number of raw files in a day. If the default is not working, look at numbers printed in the command window from step 1
to guess what a typical number of raw files per day is. 
     av = axis limits in Hz or kHz, depending on decimation factor.

I haven't really needed to mess with the other variables, read through the code if you're curious what they do.

Run this on the output from step 1 for each deployment. 


Step 3: plotDailyAvesCombinedLSGs_step3.m

Run on all deployments together, after completing steps 1 and 2 for each deployment.
This creates the final LSG.

Sean has a clear vision about how these plots should look. Changes I've made to improve appearance for some sample rates and 
data durations may negatively impact appearance for other frequencies, and could inspire outrage. You've been warned ;)