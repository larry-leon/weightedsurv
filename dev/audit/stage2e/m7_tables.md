### M7 requirements (m7_required.txt)

```
FAIL 1 (cells 7, 8, 9, 15, 19, replicates 1-500; package "estimated" against the evaluation's script: sig2_surv0, sig2_surv1, sig2_dhat at the time points, band critical value and limits, RMST variance, critical value and band limits, sig_bhat_asy and |z| for FH(0,0) and FH(0,1)): 3 of 2500 replicates failing; largest difference 0.000418 (limit 1e-8); the failing replicates are in the Cox / log-rank part only, see m7_req1_diagnosis.txt
OK 2 (cells 31, 32; "stabilized" against "ipw" from the same model, both options: sig2_dhat, band critical value, RMST variance): 0 of 4000 replicates failing; largest difference 4.72e-14 (limit 1e-8)
OK 3 (every cell; dhat, rmst_time, bhat and lr identical under "known" and "estimated"): 0 of 22500 replicates failing; 0 comparisons differ; 15 cells
OK 4 (cell 34, tied times; "estimated" FH(0,0) against the coxph residual routes, Breslow): 0 of 2000 replicates failing; largest difference of z^2 1.56e-10 (limit 1e-8); largest difference of the SE (limit 1e-6): 3.5e-10
```

Notes: 
- requirement 1, cell 7: 0 of 500 replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities 2.09e-14, of the Cox / log-rank quantities 3.64e-12; the script's alignment confirmation passed in 500 (FH(0,0)) and 500 (FH(0,1)) of 500 replicates
- requirement 1, cell 8: 0 of 500 replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities 3.62e-14, of the Cox / log-rank quantities 1.09e-12; the script's alignment confirmation passed in 499 (FH(0,0)) and 499 (FH(0,1)) of 500 replicates
- requirement 1, cell 9: 1 of 500 replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities 6.86e-14, of the Cox / log-rank quantities 0.000199; the script's alignment confirmation passed in 499 (FH(0,0)) and 499 (FH(0,1)) of 500 replicates
- requirement 1, cell 15: 2 of 500 replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities 7.33e-14, of the Cox / log-rank quantities 0.000418; the script's alignment confirmation passed in 500 (FH(0,0)) and 500 (FH(0,1)) of 500 replicates
- requirement 1, cell 19: 0 of 500 replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities 3.71e-14, of the Cox / log-rank quantities 2.53e-14; the script's alignment confirmation passed in 500 (FH(0,0)) and 500 (FH(0,1)) of 500 replicates

### Table M7 — operating characteristics with the package code (nominal 0.95; two-sided 5% log-rank test)

Rows: cell x option (analyses A, B, C in cell 40; the two bootstrap rows in cells 38 and 39, first 500 replicates). ESS ratio = mean of ess / n by arm (control / treated).
Length ratio = mean interval length under "estimated" over that under "known" (shown on the "estimated" row). `—` = not evaluated.

| k | Cell | Option | ESS ratio c / t | delta(18) | band [3, 42] | RMST(36) | log HR FH(0,0) | log HR FH(0,1) | SE/SD FH(0,0) | reject FH(0,0) | reject FH(0,1) | length ratio delta(18) | length ratio RMST(36) | length ratio log HR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | E-200 | known (500) | 0.84 / 0.70 | 0.966 | 0.962 | 0.962 | 0.970 | 0.948 | 1.120 | — | — | — | — | — |
| 7 | E-200 | estimated (500) | 0.84 / 0.70 | 0.922 | 0.934 | 0.934 | 0.934 | 0.926 | 0.948 | — | — | 0.869 | 0.828 | 0.847 |
| 8 | E-400 | known (500) | 0.82 / 0.69 | 0.972 | 0.978 | 0.976 | 0.980 | 0.962 | 1.189 | — | — | — | — | — |
| 8 | E-400 | estimated (500) | 0.82 / 0.69 | 0.934 | 0.956 | 0.934 | 0.956 | 0.952 | 1.002 | — | — | 0.863 | 0.820 | 0.843 |
| 9 | E-800 | known (500) | 0.83 / 0.69 | 0.962 | 0.988 | 0.972 | 0.980 | 0.960 | 1.125 | — | — | — | — | — |
| 9 | E-800 | estimated (500) | 0.83 / 0.69 | 0.936 | 0.960 | 0.926 | 0.944 | 0.944 | 0.952 | — | — | 0.866 | 0.821 | 0.847 |
| 15 | E-400-null | known (500) | 0.83 / 0.70 | 0.964 | 0.962 | 0.966 | 0.958 | 0.948 | 1.040 | 0.044 | 0.064 | — | — | — |
| 15 | E-400-null | estimated (500) | 0.83 / 0.70 | 0.938 | 0.950 | 0.928 | 0.940 | 0.914 | 0.879 | 0.070 | 0.092 | 0.870 | 0.826 | 0.845 |
| 19 | E-400-delay | known (500) | 0.84 / 0.70 | 0.980 | 0.984 | 0.984 | — | — | 1.153 | — | — | — | — | — |
| 19 | E-400-delay | estimated (500) | 0.84 / 0.70 | 0.944 | 0.966 | 0.954 | — | — | 0.974 | — | — | 0.865 | 0.824 | 0.845 |
| 31 | Es-400 | known | 0.83 / 0.69 | 0.970 | 0.976 | 0.981 | 0.976 | 0.952 | 1.154 | — | — | — | — | — |
| 31 | Es-400 | estimated | 0.83 / 0.69 | 0.941 | 0.949 | 0.946 | 0.944 | 0.936 | 0.976 | — | — | 0.865 | 0.821 | 0.846 |
| 32 | Es-400-null | known | 0.83 / 0.69 | 0.971 | 0.962 | 0.976 | 0.968 | 0.940 | 1.099 | 0.035 | 0.067 | — | — | — |
| 32 | Es-400-null | estimated | 0.83 / 0.69 | 0.943 | 0.939 | 0.945 | 0.943 | 0.925 | 0.930 | 0.060 | 0.085 | 0.870 | 0.826 | 0.847 |
| 33 | Es-400-null-strata | known | 0.83 / 0.69 | — | — | — | — | — | — | 0.057 | 0.076 | — | — | — |
| 33 | Es-400-null-strata | estimated | 0.83 / 0.69 | — | — | — | — | — | — | 0.065 | 0.081 | — | — | — |
| 34 | Et-400-null | known | 0.83 / 0.69 | 0.973 | 0.963 | — | 0.974 | 0.951 | 1.161 | 0.028 | 0.056 | — | — | — |
| 34 | Et-400-null | estimated | 0.83 / 0.69 | 0.948 | 0.935 | — | 0.940 | 0.935 | 0.979 | 0.060 | 0.070 | 0.867 | — | 0.844 |
| 35 | Fe-400 | known | 0.50 / 0.39 | 0.959 | 0.897 | 0.957 | 0.933 | 0.898 | 1.006 | — | — | — | — | — |
| 35 | Fe-400 | estimated | 0.50 / 0.39 | 0.916 | 0.864 | 0.916 | 0.901 | 0.880 | 0.848 | — | — | 0.850 | 0.822 | 0.843 |
| 36 | Fe-800 | known | 0.47 / 0.33 | 0.963 | 0.920 | 0.951 | 0.936 | 0.909 | 0.946 | — | — | — | — | — |
| 36 | Fe-800 | estimated | 0.47 / 0.33 | 0.925 | 0.892 | 0.903 | 0.894 | 0.892 | 0.807 | — | — | 0.855 | 0.831 | 0.854 |
| 37 | Fe-400-null | known | 0.51 / 0.37 | 0.942 | 0.865 | 0.926 | 0.914 | 0.864 | 0.897 | 0.095 | 0.143 | — | — | — |
| 37 | Fe-400-null | estimated | 0.51 / 0.37 | 0.905 | 0.818 | 0.883 | 0.871 | 0.847 | 0.760 | 0.136 | 0.164 | 0.860 | 0.830 | 0.847 |
| 38 | V-1450-null | known | 0.95 / 0.36 | 0.943 | 0.875 | 0.931 | 0.906 | — | 0.925 | 0.103 | — | — | — | — |
| 38 | V-1450-null | estimated | 0.95 / 0.36 | 0.925 | 0.840 | 0.910 | 0.889 | — | 0.861 | 0.117 | — | 0.933 | 0.916 | 0.931 |
| 38 | V-1450-null | boot (500) | — | 0.906 | — | 0.898 | 0.882 | — | 0.846 | — | — | — | — | — |
| 38 | V-1450-null | boot_pct (500) | — | 0.912 | — | 0.906 | 0.876 | — | 0.846 | — | — | — | — | — |
| 39 | V-1450-ph | known | 0.95 / 0.35 | 0.948 | 0.932 | 0.946 | 0.920 | — | 0.962 | — | — | — | — | — |
| 39 | V-1450-ph | estimated | 0.95 / 0.35 | 0.931 | 0.897 | 0.926 | 0.902 | — | 0.893 | — | — | 0.920 | 0.905 | 0.928 |
| 39 | V-1450-ph | boot (500) | — | 0.932 | — | 0.920 | 0.892 | — | 0.889 | — | — | — | — | — |
| 39 | V-1450-ph | boot_pct (500) | — | 0.932 | — | 0.922 | 0.894 | — | 0.889 | — | — | — | — | — |
| 40 | Sub-800-null | A | 0.86 / 0.95 | 0.948 | — | — | 0.942 | 0.942 | 0.985 | 0.059 | 0.065 | — | — | — |
| 40 | Sub-800-null | B | 0.86 / 0.95 | 0.948 | — | — | 0.951 | 0.941 | 0.986 | 0.053 | 0.066 | 0.957 | — | 0.927 |
| 40 | Sub-800-null | C | 0.86 / 0.95 | 0.958 | — | — | 0.966 | 0.947 | 1.064 | 0.036 | 0.059 | — | — | — |

### Table M7-B — bootstrap benchmark, the same first 500 replicates (coverage; SE/SD in brackets)

| k | Cell | Quantity | known | estimated | boot | boot_pct |
|---|---|---|---|---|---|---|
| 38 | V-1450-null | delta(18) | 0.930 (0.96) | 0.912 (0.90) | 0.906 (0.91) | 0.912 (0.91) |
| 38 | V-1450-null | R(36) | 0.926 (0.96) | 0.900 (0.88) | 0.898 (0.89) | 0.906 (0.89) |
| 38 | V-1450-null | log HR | 0.898 (0.89) | 0.884 (0.83) | 0.882 (0.85) | 0.876 (0.85) |
| 39 | V-1450-ph | delta(18) | 0.948 (1.01) | 0.932 (0.94) | 0.932 (0.95) | 0.932 (0.95) |
| 39 | V-1450-ph | R(36) | 0.936 (1.01) | 0.920 (0.91) | 0.920 (0.93) | 0.922 (0.93) |
| 39 | V-1450-ph | log HR | 0.914 (0.94) | 0.890 (0.87) | 0.892 (0.89) | 0.894 (0.89) |

### Correlation of the log-rank score with the survival difference at month 18 (FH(0,0), "estimated")

| k | Cell | mean estimated correlation (wlr_dhat_estimates) | empirical correlation over replicates |
|---|---|---|---|
| 31 | Es-400 | 0.708 | 0.722 |
| 32 | Es-400-null | 0.737 | 0.725 |

### Failures and timing

| k | Cell | replicates | ps_weights() errors | largest count of failed calls (any quantity) | largest count of non-finite intervals | seconds | cores | package SHA |
|---|---|---|---|---|---|---|---|---|
| 7 | E-200 | 500 | 0 | 0 | 1 | 13 | 12 | c195953 |
| 8 | E-400 | 500 | 0 | 0 | 0 | 23 | 12 | c195953 |
| 9 | E-800 | 500 | 0 | 0 | 0 | 60 | 12 | c195953 |
| 15 | E-400-null | 500 | 0 | 0 | 0 | 23 | 12 | c195953 |
| 19 | E-400-delay | 500 | 0 | 0 | 0 | 24 | 12 | c195953 |
| 31 | Es-400 | 2000 | 0 | 0 | 0 | 113 | 12 | c195953 |
| 32 | Es-400-null | 2000 | 0 | 0 | 0 | 108 | 12 | c195953 |
| 33 | Es-400-null-strata | 2000 | 0 | 0 | 0 | 30 | 12 | c195953 |
| 34 | Et-400-null | 2000 | 0 | 0 | 0 | 43 | 12 | c195953 |
| 35 | Fe-400 | 2000 | 0 | 0 | 0 | 73 | 12 | c195953 |
| 36 | Fe-800 | 2000 | 0 | 0 | 0 | 172 | 12 | c195953 |
| 37 | Fe-400-null | 2000 | 0 | 0 | 0 | 70 | 12 | c195953 |
| 38 | V-1450-null | 2000 | 0 | 0 | 0 | 504 | 12 | c195953 |
| 39 | V-1450-ph | 2000 | 0 | 0 | 0 | 480 | 12 | c195953 |
| 40 | Sub-800-null | 2000 | 0 | 0 | 0 | 58 | 12 | c195953 |

