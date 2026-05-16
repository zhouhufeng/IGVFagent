# Network-integration demo — synthetic signed cascade

- PKN: 6 edges (EGFR → SOS1 → RAS → RAF → MEK → ERK → MYC)
- Perturbation: EGFR=+1   Measurement: MYC=+1.2
- Status: **optimal**, objective=0.30000000000000004
- Selected: **6** edges (beta=0.05, solver=SCIP)
- Warehouse edges added: **6** (relation=`activates`, upstream=`network:demo`)

## Subnetwork
| src | sign | dst | role |
|---|---|---|---|
| EGFR | 1 | SOS1 | activates |
| SOS1 | 1 | RAS | activates |
| RAS | 1 | RAF | activates |
| RAF | 1 | MEK | activates |
| MEK | 1 | ERK | activates |
| ERK | 1 | MYC | activates |
