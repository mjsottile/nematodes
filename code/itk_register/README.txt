----------------------------
- Worm registration code   -
- Release date: 05/17/2009 -
----------------------------

Steps to use:

1) Determine the fully qualified path of the directory containing your frames.
   You do not need to provide a trailing slash.

     Example:  /Users/matt/worms/my_worm_frames

2) Create an empty directory to hold the output frames after registration
   and splitting.

     Example:  /Users/matt/worms/my_outputs

3) Modify the included input parameter file (sample.in) to reflect the
   paths, filenames, frame counts, and pixel thresholds for your frame
   sequence.  You should create one input file per movie that you are
   registering.

4) Execute the code.

     % ./wormRegister myinputfile.in

5) Wait patiently.

6) When the code completes, you will fine the output directory
   populated with a set of files.  For each input frame
   "frameXXXX.tif", you will have two files "left_frameXXXX.tif" and
   "right_frameXXXX.tif".  You will also see a file with an extension
   ".params".  This is the sequence of registration parameters used to
   transform each frame to match the reference.  These can be used to
   replay the registration, and are used when you restart at a middle
   frame to determine what the sequence of transformations up to the
   restart point were.



Please send any bugs/errors/issues to matt@cs.uoregon.edu.

