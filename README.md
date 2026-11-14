# Identification of Time-Varying Financial Networks

Dynamic network analysis of 29 oil-sensitive U.S. listed companies around the Russia-Ukraine war (24 February 2022), using system identification in MATLAB.

> **Group project** (team of three) developed for the course *System Identification and Data Analysis*, MSc in Financial Engineering, University of Padova. The report and the code were produced collaboratively.

## Goal

Understand how the web of dependencies between companies exposed to oil prices (airlines, cruise and travel, energy, logistics) changes around a geopolitical shock, and measure how oil affects each company over time.

## Data

The dataset is **not included** in this repository. The script expects an Excel file `dataset_MIMO_30.xlsx` with:

| Column | Content |
|---|---|
| 1 | `Date` (format `yyyy-MM-dd HH:mm:ss`) |
| 2 | Exogenous input `u`: daily log-returns of crude oil futures (`CL=F`) |
| 3-31 | Outputs `y`: daily log-returns of 29 listed companies |

- **Frequency:** daily, trading days only, about 1,200 observations (2021-2026), prices in USD.
- **Companies (29):** AAL, ALGT, ALK, BKNG, CCL, COP, CPA, CSX, CVX, DAL, EOG, EXPE, FDX, JBLU, LUV, MAR, MPC, NCLH, NSC, OXY, RCL, SKYW, SLB, SNCY, UAL, UNP, UPS, VLO, XOM.
- **Sectors:** airlines, cruise and travel, energy, logistics and rail.
- **Preprocessing:** log-returns r_t = ln(P_t / P_{t-1}), missing values filled, then detrending.

Market data come from third-party providers whose terms do not allow redistribution, so you need to download the prices yourself and build the file in the format above.

## Method

1. **MISO ARX models on sliding windows.** For each company *i* and each 120-day window (step 1 day), the returns of company *i* are modelled with an ARX model whose inputs are the oil return and the returns of the other 28 companies (`na = nb = 30`, `nk = 1`). Each window is normalised by its standard deviation and detrended.
2. **Regularised least squares with SS kernel.** With many parameters per node and few observations per window, the impulse responses are regularised with a stable-spline (SS) kernel estimated by `arxRegul`, which controls variance and avoids overfitting.
3. **Iterative reweighting.** A second reweighted step penalises each block of the impulse response by its norm, giving sparser and more interpretable networks.
4. **Network extraction.** The norm of each impulse-response block gives the weight of the link between two companies (Y to Y), of the oil effect on each company (U to Y), and of the self-memory of each company (AR part).
5. **Dynamic threshold.** The binary adjacency matrix is obtained with a threshold that changes in each window to keep a target density of 8.5% of the off-diagonal links.
6. **Shock date.** A network-stability metric (norm of the difference between consecutive weighted networks) is computed over time. The identified shock date is its maximum within 20 days of the real event.
7. **Oil impact and reaction time.** The mean oil impact over time is tracked, and the reaction time is the delay between the shock and the peak of oil impact within a 60-day window.
8. **Validation.** Hold-out cross-validation: the first 75% of the data for training and the last 25% for validation. The one-step-ahead fit (`compare`, zero initial conditions) is computed for each company and averaged.

## Outputs

The script produces:
- adjacency matrices of the network from 90 days before to 90 days after the shock, with zoom on the pre- and post-shock periods;
- oil impact (U to Y) and self-memory (AR) per company over time, as heatmaps and ranked bar plots by sector;
- the network-stability curve with real and identified shock dates, and the evolution of the dynamic threshold;
- the oil-impact curve with the peak after the shock and the reaction time in trading days;
- the hold-out fit for each company and its average.

Figures are included in `financial_network_identification.pdf`.

## How to run

**Requirements:** MATLAB with the System Identification Toolbox and the Statistics and Machine Learning Toolbox (`prctile`).

1. Put `dataset_MIMO_30.xlsx` in the same folder as the script.
2. Open `financial_network_identification.m` (or the Live Script) and run it section by section.
3. The shock date, the identified date and the reaction time are printed in the command window.

## Repository content

```
financial-network-identification/
├── README.md
├── financial_network_identification.m      MATLAB script
├── financial_network_identification.pdf    exported figures and results
└── .gitignore
```

## Notes and limitations

- The network shows statistical dependence between returns, not causal effects.
- Results depend on the window length (120 days), the model order (30 lags), the target density (8.5%) and the reaction window (60 days).
- The identified shock date depends on the metric chosen and on the search margin around the real event.
- The one-step-ahead fit on the validation set is modest in absolute terms, which is expected with noisy daily returns. The result shows that the regularisation prevents overfitting, not that returns are predictable.

## Authors

Group project by three students of the MSc in Financial Engineering, University of Padova. My contribution: [describe here what you did].
