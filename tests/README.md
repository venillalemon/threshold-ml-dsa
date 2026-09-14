# Mixed-security backend tests

Standalone networked tests for the GMW (offline circuit) + WRK (online circuit)
backend. Each runs a small circuit across nP local processes.

- `backend_smoke_test.cpp` — GMW COT setup, `random()`, checked `open()`.
- `gmw_eval_test.cpp` — `gmw.evaluate()` of an 8-bit adder (offline circuit path).
- `wrk_bridge_test.cpp` — `wrk_offline` (garble + open output masks offline)
  then `wrk_online` (one flight of labels + local decode) of a toy `C_post`.

Build with `-DMLDSA_TESTS=ON`, then run each `PARTY` in 1..nP with a shared
`EMP_PORT`, e.g. `EMP_PORT=17720 ./build/wrk_bridge_test3 1` for parties 1..3.
