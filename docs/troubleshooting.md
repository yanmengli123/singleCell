# Troubleshooting Guide

## Environment Issues

### torch.cuda Warning
```
WARNING: The pynvml package is deprecated. Please install nvidia-ml-py instead.
```
- **Severity**: Low (warning only)
- **Cause**: `pynvml` is outdated, `nvidia-ml-py` is recommended
- **Solution**: This is non-critical and does not affect functionality

### TensorFlow Warnings
```
WARNING: tf.losses.sparse_softmax_cross_entropy is deprecated
```
- **Severity**: Low (warning only)
- **Cause**: TensorFlow/Keras version mismatch
- **Solution**: Non-critical for our pipelines

### oneDNN Warning
```
I tensorflow/core/util/port.cc:113] oneDNN custom operations are on.
```
- **Severity**: Low (informational)
- **Cause**: TensorFlow oneDNN optimizations enabled
- **Solution**: Normal behavior, can be disabled with `TF_ENABLE_ONEDNN_OPTS=0`

## Data Issues

### GBK Encoding Error
```
UnicodeEncodeError: 'gbk' codec can't encode character '\U0001f449'
```
- **Severity**: High (blocks output)
- **Cause**: Windows console uses GBK encoding, emoji characters can't display
- **Solution**: Remove emoji/Chinese characters from print statements, use UTF-8

### FutureWarning for __version__
```
FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')`
```
- **Severity**: Low (deprecation warning)
- **Cause**: Package version access pattern changed
- **Solution**: Update to use `importlib.metadata.version()` in future

## Performance Issues

### Slow Scanpy Import
- **Severity**: Low (once per session)
- **Cause**: scanpy loads many dependencies on first import
- **Solution**: Normal for first import, subsequent imports are faster

### Sparse Matrix Densification Warning
```
UserWarning: zero-centering a sparse array/matrix densifies it.
```
- **Severity**: Low (informational)
- **Cause**: PCA/UVM with sparse data requires densification
- **Solution**: Expected behavior for sparse preprocessing

## Leiden Algorithm Warning
```
FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.
```
- **Severity**: Medium (will change in future)
- **Cause**: Leiden library change
- **Solution**: Pass `flavor="igraph"` and `n_iterations=2` to prepare for future

## Common Fixes

1. **UTF-8 Output**: Always use UTF-8 compatible strings in output
2. **File Paths**: Use `pathlib.Path` for cross-platform compatibility
3. **Random State**: Always set `random_state=0` for reproducibility
4. **Memory**: For large datasets, process in batches if needed