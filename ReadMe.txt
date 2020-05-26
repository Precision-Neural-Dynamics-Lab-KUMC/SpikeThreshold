
extractSpikesFromBroadband takes continuously collected neural data in a .ns5 file 
and calculates a threshold based on the signal median and outputs spike snippets in a .nex file for further spike sorting

Required tools for extractSpikesFromBroadband:
1) NPMK for reading .ns5 files
Available here:
https://github.com/BlackrockMicrosystems/NPMK/tree/master/NPMK


2) HowToReadAndWriteNexAndNex5FilesInMatlab from NeuroExplorer
Available here:  
Code to Read and Write NeuroExplorer Data Files, Matlab code to read and write .nex files
https://www.neuroexplorer.com/downloadspage/

*NOTE, the writeNexFile function does not create files that read into Plexon Offline Sorter the same as the native export/import .NEX in Plexon 
A modified writeNexFile.m is included with WaveLimit that is recommended for use with Plexon Offline Sorter.  
Simply replace writeNexFile.m in the HowToReadAndWriteNexAndNex5FilesInMatlab directory. 

To run, NPMK and HowToReadAndWriteNexAndNex5FilesInMatlab must be on your Matlab path.
You can either add it permanantly or include it using the example calls in WaveLimit_example_call_script.m