im_start = 1;
im_end = 7200;
im_ref = 1;

postmask = [1 160 1 640; 1 220 270 282; 1 220 298 305];
premask = [1 480 1 190;1 480 390 640; 390 480 1 640];

fname_base = '/Volumes/My Passport/worms/fall09/fall09 ';

[r,c,h]=process(fname_base,'.png', im_start, im_end, im_ref, 165, premask, postmask, 1);
