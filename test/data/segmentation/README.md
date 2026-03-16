# Segmentation

For segmentation, we provide pre-compiled binaries of RawHash1 and RawHash2 which only perform the segmentation step, using the Scrappie (R9.4) and Scrappie (R10.4.1) segmentation algorithms respectively.

rawhash_segment_signal runs just as RawHash would, but it prints out the events in the following format:
```
read_id1 event1 event2 event3 ... eventn 
read_id2 event1 event2 event3 ... eventn 
```

To get some examples on how to use it, check out the event generation scripts in /CERN/test/data/generate_events