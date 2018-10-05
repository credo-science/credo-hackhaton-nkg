# credo-hackathon-nkg

Extract and simple classify cosmic-ray hits from FITS files with dark frames. Supports removing hot pixels
based on next frame when is necessary.

The algorithm:
1. (optional) remove hot pixels based on next frame by `ccdproc` library.
2. Finding hits by `source_properties` function from `photutils` library.
3. Classify hits to `dots`, `comets` and `worms` using `ellipticity` parameter from `source_properties` function:
   *  when `ellipticity <= 0.2` then classified as `dot`,
   *  when `ellipticity >= 0.6` then classified as `comet`,
   *  otherwise `ellipticity <= 0.2` was classified as `worm`.
4. Detect group of hits. Hits are in groups where distance of pair of hits is less than 30px.
5. Counting all groups by *number of groups formed in a graph of friends* algorithm.
6. Print stats for each FITS file contains count of hits and area of hits by each class
and count of groups.

## How to prepare to run

```bash
virtualenv -p python3 env
source env/bin/activate
pip install -r requirements.txt
```

## How to run

1. Unpack FITS with dark frames with hot pixels or dark frames with hot pixels
removed to `image_dir` directory.
2. In file `filter_ccd.py` change:
   * `image_dir` - path to directory where you unpack FITS files,
   * `output_dir` - path where product of algorithm will be stored,
   * `files` - array of array of file name with frames with hot pixels,
   one array is the series of on observation, hot pixels will be removed
   based on next file, finding of hits is doing on 1 to N-1 files from the array,
   * `darks` - array of files with hit pixels removed.
3. Run: 
```bash
python filter_ccd.py
```

## How to interpret result

### Values printed to `stdout`
The `filter_ccd.py` print *gnuplot's* friendly `dat` file format output to `stdout`.

Columns:
1. FITS file name
2. count of hits classified as dot
3. area of hits classified as dot
4. count of hits classified as comet
5. area of hits classified as comet
6. count of hits classified as worm
7. area of hits classified as worm
8. count of all hits
9. area of all hits
10. count of hits in groups
11. area of hits in groups
12. count of all groups
13. exposure time from FITS metadata
14. date of start of exposure from FITS metadata
15. and other metadata from FITS configured in `metadata` const in `filter_ccd.py` 

### Files in `output_dir` 
In `output_dir` you get:
* `cleared-{file_from_files}.fits` - FITS file with removed hot pixels
* `cleared-{file_from_files_or_darks}.png` - PNG export with pixels bright mapped from
`[A,B]` to `[0,255]` by linear function and `A` is percentile 25%, `B` is percentile 99.999%,
* `detections/{file_from_files_or_darks}/` - result of detection:
   *  `./{a}-{b}-{c}-{d}-{e}-{f}.fits` - cut hit, where
   `a` - area of hit in pixels,
   `b` - XY of hit,
   `c` - class of hit,
   `d` - ellipticity,
   `e` - elongation,
   `f` - eccentricity,
   * `./{a}-{b}-{c}-{d}-{e}-{f}.png` - PNG export of above from `cleared-{file_from_files_or_darks}.png`,
* `detections/{file_from_files_or_darks}/classes/{a}-{b}-{c}-{d}-{e}-{f}.png` - marked area by *find_hits* function.
