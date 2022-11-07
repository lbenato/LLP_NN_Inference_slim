# Application of TensorFlow interface in CMSSW to LLP NN Inference

A small example is included in folder ```nn_inference```: a keras model and a small tree with 100 events are provided.

### Setup

Install CMSSW release and cmsml package [https://cms-ml.github.io/documentation/inference/tensorflow2.html#saving-your-model](https://cms-ml.github.io/documentation/inference/tensorflow2.html#saving-your-model)

Note: SL7 is required

```bash
cmsrel CMSSW_12_5_2
cd CMSSW_11_5_2/src
cmsenv
git cms-init
scram b -j 32
cd $CMSSW_BASE/src
pip install --upgrade --user git+https://github.com/cms-ml/cmsml
```

Clone and compile this repo

```bash
mkdir NNInferenceCMSSW_slim
cd NNInferenceCMSSW_slim
git clone https://github.com/lbenato/LLP_NN_Inference_slim.git
scram b -j 32
```

### Run NN inference
```bash
cd $CMSSW_BASE/src
../bin/slc7_amd64_gcc820/tf_test
```