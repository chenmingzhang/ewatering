# ewatering
modelling scripts for e-watering project



```
grep -rnw '/mnt/d/projects/flopy/examples' -e 'ModflowRch'
chenming@modflow$ grep -rnw '/mnt/d/projects/flopy/examples' -e 'ModflowRch'
/mnt/d/projects/flopy/examples/groundwater_paper/Notebooks/example_1.ipynb:674:    "fpm.ModflowRch(model, rech=0.001)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_swi2package_ex4-checkpoint.ipynb:342:    "rch = flopy.modflow.ModflowRch(ml, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_swi2package_ex4-checkpoint.ipynb:440:    "rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_WatertableRecharge_example-checkpoint.ipynb:261:    "rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=1)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:113:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:368:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:711:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_modpath7_structured_example.ipynb:205:    "flopy.modflow.ModflowRch(m, ipakcb=iu_cbc, rech=rch)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_MT3DMS_examples.ipynb:1709:    "    rch = flopy.modflow.ModflowRch(mf, rech=rech)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_MT3DMS_examples.ipynb:2207:    "    rch = flopy.modflow.ModflowRch(mf, rech=1.14e-3)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_multi-component_SSM.ipynb:102:    "rch = flopy.modflow.ModflowRch(mf)"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_PEST.ipynb:529:    "rch = flopy.modflow.ModflowRch(m, rech={0: 0.001, 2: 0.003})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_swi2package_ex4.ipynb:342:    "rch = flopy.modflow.ModflowRch(ml, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_swi2package_ex4.ipynb:440:    "rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_WatertableRecharge_example.ipynb:261:    "rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=1)\n",
/mnt/d/projects/flopy/examples/scripts/flopy_swi2_ex4.py:210:        rch = flopy.modflow.ModflowRch(ml, rech=rch_data)
/mnt/d/projects/flopy/examples/scripts/flopy_swi2_ex4.py:267:        rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)chenming@modflow$ grep -rnw '/mnt/d/projects/flopy/examples' -e 'ModflowRch'
/mnt/d/projects/flopy/examples/groundwater_paper/Notebooks/example_1.ipynb:674:    "fpm.ModflowRch(model, rech=0.001)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_swi2package_ex4-checkpoint.ipynb:342:    "rch = flopy.modflow.ModflowRch(ml, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_swi2package_ex4-checkpoint.ipynb:440:    "rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/.ipynb_checkpoints/flopy3_WatertableRecharge_example-checkpoint.ipynb:261:    "rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=1)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:113:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:368:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_external_file_handling.ipynb:711:    "rch = flopy.modflow.ModflowRch(ml, rech={0: warmup_recharge, 1: fname})\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_modpath7_structured_example.ipynb:205:    "flopy.modflow.ModflowRch(m, ipakcb=iu_cbc, rech=rch)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_MT3DMS_examples.ipynb:1709:    "    rch = flopy.modflow.ModflowRch(mf, rech=rech)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_MT3DMS_examples.ipynb:2207:    "    rch = flopy.modflow.ModflowRch(mf, rech=1.14e-3)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_multi-component_SSM.ipynb:102:    "rch = flopy.modflow.ModflowRch(mf)"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_PEST.ipynb:529:    "rch = flopy.modflow.ModflowRch(m, rech={0: 0.001, 2: 0.003})"
/mnt/d/projects/flopy/examples/Notebooks/flopy3_swi2package_ex4.ipynb:342:    "rch = flopy.modflow.ModflowRch(ml, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_swi2package_ex4.ipynb:440:    "rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)\n",
/mnt/d/projects/flopy/examples/Notebooks/flopy3_WatertableRecharge_example.ipynb:261:    "rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=1)\n",
/mnt/d/projects/flopy/examples/scripts/flopy_swi2_ex4.py:210:        rch = flopy.modflow.ModflowRch(ml, rech=rch_data)
/mnt/d/projects/flopy/examples/scripts/flopy_swi2_ex4.py:267:        rch = flopy.modflow.ModflowRch(ml2, rech=rch_data)

```


fhb package


```
chenming@modflow$ grep -rnw '/mnt/d/projects/flopy/examples' -e 'ModflowFhb'
chenming@modflow$

```


```
enming@modflow$ grep -rnw '/mnt/d/projects/flopy/examples' -e 'ModflowLak'
/mnt/d/projects/flopy/examples/Notebooks/flopy3_mt3d-usgs_example_with_sft_lkt_uzt.ipynb:529:    "lak = flopy.modflow.ModflowLak(\n",
```
