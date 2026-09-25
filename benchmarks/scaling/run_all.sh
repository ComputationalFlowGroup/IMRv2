#!/bin/bash

# Copyright (C) The IMR developers (see README.md)
# SPDX-License-Identifier: GPL-3.0-or-later

sbatch Rresolved.sh
sbatch NMt_resolution.sh
sbatch Nt_resolution.sh
