# Overview
`TrigSync` can be used to find a set of trigger pulses from one electrophysiological recording inside another. The compatible file formats are Blackrock (NSx, NEV) and the clinical European Data Format (EDF).

# Requirements
This package has been tested on Windows 10 with MATLAB 2018b, though generally versions 2017b and later should be sufficient. You will also need to have the following libraries accessible in the path:
* [NPMK](https://github.com/BlackrockMicrosystems/NPMK) (tested with v4.4.2.0)
* [FieldTrip](http://www.fieldtriptoolbox.org/) (tested with v20160929)

# Installation
After extracting the package to a directory of your choosing, add the directory and its contents by running `addpath(genpath(*your-path-here*))`.

# :construction: Operation
The key function in this package is `locateCodedTrigger`. The function can be used to find a set of coded triggers in one *target* recording, within another *origin* recording. A coded trigger is a set of trigger pulses whose value encodes information about its relative location in the recording.

## locateCodedTrigger
This function returns a structure that holds two items:
1. the first sample of the first trigger pulse found in the target recording, and
2. the corresponding sample of that same trigger pulse within the origin recording.

```MATLAB
t0 = locateCodedTrigger(targetpath,originpath,'presentation','targettrigchan',101,'targetoffsettime',offsettime,'validate',true);
```

Here, `locateCodedTrigger` is looking for `presentation`-type triggers (these encode hours/minutes/seconds using delays between sets of pulses). The parameter `targetoffsettime` (double) specifies in seconds a temporal offset that can be used to look for trigger pulses that are somewhere in the middle of the target recording.

Setting `validate=true` produces a plot that overlays the trigger pulses from both target and origin files, such that they may be inspected.

<p align="center">
    <img src="https://github.com/Center-For-Neurotechnology/trigsync/raw/master/images/validation.png"/><br>
    <b>Fig. 1.</b> Overlays of target NSx (blue) and origin EDF (orange) trigger pulses.
</p>

# Note
This version of `TrigSync` is fresh, and likely bug-ridden! Please keep this in mind if you want to use it for your analysis. Any feedback is very welcome.
