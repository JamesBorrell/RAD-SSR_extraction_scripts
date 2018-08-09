## PIPELINE PLAN
## 
### 1. Define the repeat motifs you would like to search for, simple bialleleic repeats are a good place to start.
### 2. Search for RAD tags that contain these motifs
### 3. Filter the selected tags to retain the best candidates. For example, those present in many individuals and without repeats at the beginning and end of the read etc.
### 4. Generate 'primers' for those RAD tags, in the form of the first 15-25 bases, to allow us to search for repeat motifs in the whole read library.
### 5. Search for selected reads across all individuals.
### 6. Run a script to convert the repeat motifs in raw reads into factors, and to call heterozygotes or homozygotes depending on the number of matching reads.
### 7. Convert this into desired file format e.g. .stru
### 8. Apply any biases, if you wish, for example selecting highly variable loci in the same manner that PCR-SSRs are selected.
