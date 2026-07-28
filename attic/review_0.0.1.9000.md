# Code review: libcmaesr (whole package)

Review of the full package (R code, C/C++ wrapper, tests, build files) at commit `75bae0e`.
The two most severe findings were verified empirically with reproduction scripts, not just by reading the code.

## Summary

The package is small and the overall architecture is sound — the batch-eval delegation trick (`evalf` → cache → `strat.eval`) is genuinely clever, the exception-based error bridging between R and C++ is mostly done right, and test coverage of the algorithm matrix is thorough. But it is shipping **two confirmed correctness bugs** (silent wrong results on nested calls; R-API `longjmp` through the C++ optimizer stack on integer returns), plus a known-broken name lookup with a `FIXME: bug` comment and the correct implementation sitting commented out directly above it. Not CRAN-ready in this state.

## Critical issues (blocking)

1. **[FIXED]** **Nested `cmaes()` calls silently optimize the wrong function.** `src/cmaes_wrap.cpp:21` (`static SEXP G_OBJ`), set at `cmaes_wrap.cpp:269`, never saved/restored. If the objective itself calls `cmaes()` (nested optimization — entirely plausible for bilevel problems or an objective that calibrates a sub-model), the inner call clobbers `G_OBJ`, and after it returns the outer optimizer keeps evaluating the **inner** objective. Reproduced: outer run on `sum(x^2)` returned `x = (0.8998, 0.9002)`, `y = 8.5e-08` — the inner objective's optimum and value. No error, no warning, just wrong answers. `G_EVAL_CACHE` and `G_IN_BATCH` have the same problem in batch mode. Either make the state re-entrant (save/restore on entry/exit, RAII guard) or detect re-entry and throw.

2. **[FIXED]** **Batch mode: integer/logical objective returns cause `Rf_error` `longjmp` through live C++ frames.** `check_numvec` (`cmaes_wrap.cpp:34`) uses `Rf_isNumeric`, which accepts `INTSXP` and `LGLSXP`. Then `cmaes_wrap.cpp:106` calls `REAL(s_y)[r]`, which raises an R error on a non-`REALSXP`. Reproduced: `fn = function(x) rep(1L, nrow(x))` → `REAL() can only be applied to a 'numeric', not a 'integer'`. That error is a `longjmp` from inside `strat.optimize()` — it skips every C++ destructor on the stack (strategy, `std::function`s, Eigen buffers) and leaves `G_IN_BATCH` stuck at `true`, which then poisons a later run's cache-miss path with a spurious "pointer mismatch" error. It is also inconsistent: single mode happily accepts `1L` via `Rf_asReal` (verified). The entire `libcmaesr_error` machinery exists precisely to prevent longjmp-through-C++, and this hole defeats it. Fix: check `TYPEOF(s_y) == REALSXP` in `check_numvec` (or coerce with `Rf_coerceVector` under protection) so the failure goes through the exception path.

3. **[FIXED]** **`RC_find_name` is a known-broken prefix match, with the fix commented out above it.** `src/rc_helpers.c:54` — `strncmp(stored, name, strlen(name))` matches any stored name that merely *starts with* the query. Line 23: `// FIXME: bug: we only do a prefix check here.....`, followed by 25 lines of a correct exact-match implementation in a comment block. Today this works only by ordering luck: `"tpa"` is found before `"tpa_dsigma"` because of its position in the `cmaes_control` list. Reorder that list, or add any field sharing a prefix, and every lookup routed through `RC_list_get_el_by_name` silently reads the wrong control parameter — no error, just a misconfigured optimizer. This is load-bearing for every single control parameter. Swap in the exact-match version and delete the comment block.

4. **[FIXED]** **`G_EVAL_CACHE` holds dangling pointers across calls.** `cmaes_wrap.cpp:104` clears the cache at the start of each batch `evalf`, but nothing clears it when a run ends. The keys are heap addresses inside a destroyed `phenocands`. In a subsequent batch run, any `cachedFF` invocation *before* the first `evalf` (elitism reinjection, uncertainty-handling re-evals) does a lookup in a cache full of dangling pointers — and if the allocator reused an address, it false-hits and silently returns a stale fitness from the previous optimization. Low probability, unbounded consequence, one-line fix: clear the cache (and reset `G_IN_BATCH`) at the top of `c_cmaes_wrap`.

## Required changes

5. **[FIXED]** **Delete the dead — and wrong — transpose loop.** `cmaes_wrap.cpp:90-92` memcpys `phenocands.row(i).data()` as if a row of a column-major Eigen matrix were contiguous. It isn't; this copies interleaved garbage. It happens to be harmless only because lines 95-98 immediately overwrite the entire buffer with the correct `X = phenocands.transpose()`. The same copy is written twice, and the first one is incorrect. Remove the loop and its comment.

6.  **[FIXED]** **`call_obj_with_error_handling_PROTECT` ignores its `s_obj` parameter** (`cmaes_wrap.cpp:44-46`) and uses `G_OBJ` instead. The batch path at line 101 carefully passes the closure-captured `s_obj`, believing it matters. It doesn't. This misleading signature is exactly how bug 1 survives review. Use the parameter, drop the global from this function.

7. **[FIXED]** **The x0 bounds error message lies.** `R/cmaes.R:259-263`: the check is inclusive (`x0 >= lower`, `x0 <= upper` — deliberately, per commit "allow x0 on lower and upper"), but the message says "'x0' must be **strictly** between 'lower' and 'upper'!" and `test_single.R:284` asserts on the wrong wording. Also inconsistent: `mlr3misc::stopf` with an explicit `::` here, plain `stop()` everywhere else, violating the project style rule.

8. **[FIXED]** **`cmaes()` trusts control-object internals it never validates.** `assert_class(control, "cmaes_control")` checks a class attribute, nothing more. `ctrl$algo = 42` reaches `STRING_ELT` on a `REALSXP` (`rc_helpers.c:63`) → R-internal error longjmping through C++, same UB class as item 2. The tests themselves mutate control objects after construction (`ctrl$quiet = TRUE`), so post-hoc mutation is an anticipated pattern. Re-assert the fields consumed, or type-check on the C side.

9. **[FIXED]** **Commented-out code strewn across the codebase**: `rc_helpers.c:23-49` (27 lines), the `fevals` remnants at `cmaes_wrap.cpp:296` and `:306`, `test_single.R:4`, `test_batch.R:2`, and the roxygen `@return` entry for `fevals` at `R/cmaes.R:233-234` demoted to plain `#` comments. The `fevals` situation (upstream libcmaes miscount, per `cpp-tests/test_fevals_mismatch.cpp`) deserves a tracked issue referenced in one place, not five zombie snippets.

10. **`NA`/`NaN`/`Inf` objective returns pass unchecked.** `check_numvec` validates type and length only. An objective returning `NA_real_` feeds NaN straight into the covariance update; results silently degrade. Reject non-finite values at the boundary with a clear error naming the offending candidate, or explicitly document that non-finite returns are the user's problem.

11. **[FIXED]** **Dead helper code in a CRAN package.** `RC_tryeval_PROTECT`, `RC_intvec_*`, `RC_dblvec_copy_to_SEXP`, all four `RC_df_*` functions, `RC_r6_set_member`, `RC_set_class`, `RC_intscalar/dblscalar_create`, `RC_list_create_emptynames_PROTECT` — none are called. This is presumably a personal helper library, but for CRAN it is unreviewed attack surface and dead weight; the `RC_find_name` bug (item 3) lives in exactly this junk drawer. Trim to what's used.

12. **[FIXED]**  **Exception handling loses information and has a longjmp-over-destructor.** `c_cmaes_wrap`'s `catch (...)` (`cmaes_wrap.cpp:311`) reports "unknown error" for any `std::exception` (Eigen assert, `bad_alloc`, libcmaes-internal throw) — add `catch (const std::exception &e)` with `e.what()`. And `Rf_error` in `dispatch_with_batch_eval`'s default case (`cmaes_wrap.cpp:206`) longjmps over the live `std::function cachedFF` destructor; throw `libcmaesr_error` instead — the top-level handler exists for exactly this.

## Suggestions

13. **[FIXED]** `.Call("c_cmaes_wrap", ..., PACKAGE = "libcmaesr")` uses string lookup despite `.registration = TRUE`; use the registered native symbol object.
14. **[FIXED]** `print.cmaes_control` (`R/cmaes.R:153-158`): the variable `cc_not_default` is a lie — it filters `NULL`/`NA`, so `maximize = FALSE`, `algo = "acmaes"`, `max_fevals = 100`, `quiet = TRUE` (all defaults) always print. The inner `Filter` lambda also shadows the method's `x` argument.
15. **[FIXED]** The `seed == NA_INTEGER → 0` mapping at `cmaes_wrap.cpp:220` is dead — `cmaes()` always replaces `NA` with a random seed before the `.Call`. Also `DEBUG_PRINT("... dim: %d", dim)` at `:222` passes `R_xlen_t` to `%d` — UB the day someone flips `DEBUG_ENABLED` to 1.
16. Tests: file naming violates the project's own convention (`R/cmaes.R` → `test_cmaes.R`, not `test_single.R`/`test_batch.R`); `res_names` and the logged-sphere helper are copy-pasted between the two files (belongs in `helper.R`); `eval_log <<- rbind(eval_log, row)` is O(n²) per-eval `rbind` in the hot loop of the biggest test.
17. `cleanup` removes `src/libcmaes/*.o` but not `src/*.o` (`cmaes_wrap.o`, `init.o`, `rc_helpers.o` are sitting in `src/` right now).
18. **[FIXED]** Batch docs (`R/cmaes.R:201`) promise "a numeric matrix with `lambda` rows" — false under IPOP/BIPOP restarts where the population doubles; the test helper special-cases nine algorithms to work around the doc being wrong. Say "one row per candidate; the number of rows can vary."
19. A user interrupt landing during the R objective evaluation is swallowed by `R_tryEval` and reported as "objective evaluation failed!" rather than an interrupt.

Remaining items are minor (`cmaes_algos` living in `zzz.R` against the naming convention; `test.R`/`todo.txt` tracked in git; `max_fevals = 100` is a startlingly small default budget, though documented).

## What's genuinely good

The `RC_tryeval_nothrow` + typed-exception + single-`Rf_error`-at-the-boundary design is the correct architecture for R/C++ interop, the pointer-cache batch delegation preserves all internal libcmaes behavior instead of reimplementing it, `set_mt_feval(false)` shows the R API threading constraint was actually thought about, and the algorithm × dimension × lambda test matrix with bound-verification of every logged evaluation is more rigorous than most CRAN optimization packages. That is why the globals and the `FIXME` bug are so jarring — the hard parts are done well and the easy parts were left broken.

## Verdict

**Request changes.** Items 1-4 are confirmed or trivially-triggered correctness bugs; 1, 3, and 4 produce *silently wrong optimization results*, which for a numerical package is the worst failure class there is. All four have small, local fixes.
