# HPCC NS-3 simulator
This is an NS-3 simulator for Intelligent Computing Wide-Area RDMA project. 

It is based on NS-3 version 3.17.

## Quick Start

### Build
使用前需先配置 python2 环境并激活
```
conda create -n hpcc python=2.7
conda activate hpcc
```

`./waf configure`

Please note if gcc version > 5, compilation will fail due to some ns3 code style.  If this what you encounter, please use:

`CC='gcc-5' CXX='g++-5' ./waf configure`

### Experiment config
Please see `mix/config.txt` for example. 

`mix/config_doc.txt` is a explanation of the example (texts in {..} are explanations).

`mix/fat.txt` is the topology used in HPCC paper's evaluation, and also in PINT paper's HPCC-PINT evaluation.

### Run
The direct command to run is:
`./waf --run 'scratch/third mix/config.txt'`

./waf --run 'scratch/third mix/config_ls.txt'

We provide a `run.py` for automatically *generating config* and *running experiment*. Please `python run.py -h` for usage.
Example usage:
`python run.py --cc hp --trace flow --bw 100 --topo topology --hpai 50`

To run HPCC-PINT, try:
`python run.py --cc hpccPint --trace flow --bw 100 --topo topology --hpai 50 --pint_log_base 1.05 --pint_prob 1`

## Files added/edited based on NS3
The major ones are listed here. There could be some files not listed here that are not important or not related to core logic.

`point-to-point/model/qbb-net-device.cc/h`: the net-device RDMA

`point-to-point/model/pause-header.cc/h`: the header of PFC packet

`point-to-point/model/cn-header.cc/h`: the header of CNP

`point-to-point/model/pint.cc/h`: the PINT encoding/decoding algorithms

`point-to-point/model/qbb-header.cc/h`: the header of ACK

`point-to-point/model/qbb-channel.cc/h`: the channel of qbb-net-device

`point-to-point/model/qbb-remote-channel.cc/h`

`point-to-point/model/rdma-driver.cc/h`: layer of assigning qp and manage multiple NICs

`point-to-point/model/rdma-queue-pair.cc/h`: queue pair

`point-to-point/model/rdma-hw.cc/h`: the core logic of congestion control

`point-to-point/model/switch-node.cc/h`: the node class for switch

`point-to-point/model/switch-mmu.cc/h`: the mmu module of switch

`network/utils/broadcom-egress-queue.cc/h`: the multi-queue implementation of a switch port

`network/utils/custom-header.cc/h`: a customized header class for speeding up header parsing

`network/utils/int-header.cc/h`: the header of INT

`applications/model/rdma-client.cc/h`: the application of generating RDMA traffic

