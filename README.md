# ChromaFold

## Training
- Train model without co-accessibility component
```
python script/train_bulkOnly.py --data-path /data/leslie/gaov/project/origami_publication/data/processed_input/ -ct gm12878_hg38 umb_endo imr90
```
- Train model with co-accessibility
```
python script/train.py --data-path /data/leslie/gaov/project/origami_publication/data/processed_input/ -ct gm12878_hg38 umb_endo imr90 --mod-name bothInput
```
- Train deterministic model for full reproducibility
```
python script/train_bulkOnly.py --data-path /data/leslie/gaov/project/origami_publication/data/processed_input/ -ct gm12878_hg38 umb_endo imr90 --deterministic --mod-name deterministic
```

## Inference
