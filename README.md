# WHR Optimization Framework

[![DOI](https://zenodo.org/badge/1250919465.svg)](https://doi.org/10.5281/zenodo.20756923)

A computational framework for optimization of waste hear recovery technology that uses packed-bed direct-contact condensers.

## Overview

This framework compiles and synthesizes experimental and analytical data on heat and water recovery from waste flue gases using randomly packed bed direct-contact condensers. This synthesis automatically detects gaps in existing data to guide future experimental or numerical studies. Ultimately, the processed data serves as a training dataset for a surrogate neural network model designed for engineering and optimizing waste heat recovery (WHR) units.

The framework is at the development stage now. Currently, only the hydraulics segment is available.

<p align="center">
  <a href="https://github.com/nikitin-pro/whr-optimization-framework/wiki">
    <img alt="Documentation" src="https://img.shields.io/badge/Documentation-green">
  </a>
</p>

---

## Roadmap

- [ ] Pressure drop in randomly packed beds:
  
  - [x] Cross-validation of existing correlations.
  
  - [x] Detection of gaps in existing experimental data.
  
  - [ ] Surrogate neural network model.

- [ ] Heat-and-mass transfer in randomly packed beds:
  
  - [ ] Cross-validation of existing correlations.
  
  - [ ] Detection of gaps in existing experimental data.
  
  - [ ] Surrogate neural network model.

## Contributing

Currently contribution is limited to pressure drop models or raw data. See [CONTRIBUTING.md](CONTRIBUTING.md) for the model implementation template and registration steps. Submit new models via the [issue template](.github/ISSUE_TEMPLATE/new_model_data.yml).

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Citation

If you use this framework in your research, please cite as follows:

> Nikitin, M. (2026). WHR Optimization Framework (Version milestone). Zenodo. DOI: [10.5281/zenodo.20756923](https://doi.org/10.5281/zenodo.20756923)
