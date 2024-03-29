# Performance Analysis of OpenMPI Collective Operations on the ORFEO Cluster

## Introduction
This folder contains the final report and associated materials for the High Performance Computing 2023/2024 course exercise, focusing on the performance analysis of OpenMPI collective operations, specifically broadcast and gather, on the ORFEO cluster.

## Repository Contents
- `CORONICA_report.pdf`: The final report detailing the methodology, analysis, and findings of the performance evaluation.
- `execute_benchmarks.sh`: A useful script to configure the benchmark suite (if not already), compile and run the tested benchmarks.
- In the `graphs` folder some additional and interesting graphs can be found. The errors are calculated considering the MIN and MAX latency.

## Summary
The report presents a comprehensive evaluation of OpenMPI's broadcast and gather operations across different algorithms, leveraging the OSU Micro-Benchmarks suite. Key findings indicate significant performance variations between algorithms, emphasizing the importance of algorithm selection based on message size and network architecture.

## Methodology Overview
Performance benchmarks were conducted on the ORFEO cluster using two THIN nodes, focusing on comparing the default OpenMPI implementation with selected algorithms for broadcast and gather operations. The analysis utilized the OSU Micro-Benchmarks suite, version 7.3.
