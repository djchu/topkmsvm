#ifndef SDCA_SOLVER_EVAL_CORE_H
#define SDCA_SOLVER_EVAL_CORE_H

#include <numeric>

#include "sdca/math/blas.h"
#include "sdca/solver/data.h"

namespace sdca {

template <typename Result,
          typename Input,
          typename Output,
          template <typename, typename> class Evaluation>
inline Evaluation<Result, Output>&
eval_begin(
    dataset<Input, Output, Evaluation<Result, Output>>& d
  ) {
  // Create a new eval using the default c'tor
  d.evals.resize(d.evals.size() + 1);
  return d.evals.back();
}


template <typename Result,
          typename Input,
          template <typename, typename> class Evaluation>
inline Evaluation<Result, multiclass_output>&
eval_begin(
    dataset<Input, multiclass_output,
            Evaluation<Result, multiclass_output>>& d
  ) {
  // Create a new eval using the default c'tor
  d.evals.resize(d.evals.size() + 1);
  Evaluation<Result, multiclass_output>& eval = d.evals.back();

  // Allocate storage for the multiclass performance metric
  eval.accuracy.resize(d.out.num_classes);
  return eval;
}


template <typename Int,
          typename Result,
          typename Context>
inline void
eval_end(
    const Int num_classes,
    const Int num_examples,
    const Context& ctx,
    eval_train<Result, multiclass_output>& e
  ) {
  // Set epoch and timings
  e.epoch = ctx.epoch;
  e.solve_time_cpu = ctx.solve_time.cpu.elapsed_now();
  e.solve_time_wall = ctx.solve_time.wall.elapsed_now();
  e.eval_time_cpu = ctx.eval_time.cpu.elapsed_now();
  e.eval_time_wall = ctx.eval_time.wall.elapsed_now();

  // Compute the overall primal/dual objectives and their individual terms
  ctx.objective.update_all(
    e.primal, e.dual, e.primal_loss, e.dual_loss,
    e.primal_regularizer, e.dual_regularizer);

  // Top-k accuracies for all k
  std::partial_sum(e.accuracy.begin(), e.accuracy.end(), e.accuracy.begin());
  Result coeff(1 / static_cast<Result>(num_examples));
  sdca_blas_scal(static_cast<blas_int>(num_classes), coeff, &e.accuracy[0]);
}


template <typename Int,
          typename Result,
          typename Context>
inline void
eval_end(
    const Int,
    const Int num_examples,
    const Context& ctx,
    eval_train<Result, multilabel_output>& e
  ) {
  // Set epoch and timings
  e.epoch = ctx.epoch;
  e.solve_time_cpu = ctx.solve_time.cpu.elapsed_now();
  e.solve_time_wall = ctx.solve_time.wall.elapsed_now();
  e.eval_time_cpu = ctx.eval_time.cpu.elapsed_now();
  e.eval_time_wall = ctx.eval_time.wall.elapsed_now();

  // Compute the overall primal/dual objectives and their individual terms
  ctx.objective.update_all(
    e.primal, e.dual, e.primal_loss, e.dual_loss,
    e.primal_regularizer, e.dual_regularizer);

  e.rank_loss /= static_cast<Result>(num_examples);
}


template <typename Int,
          typename Result,
          typename Context>
inline void
eval_end(
    const Int num_classes,
    const Int num_examples,
    const Context& ctx,
    eval_test<Result, multiclass_output>& e
  ) {
  // Compute the overall primal/dual objectives and their individual terms
  ctx.objective.update_primal_loss(e.primal_loss);

  // Top-k accuracies for all k
  std::partial_sum(e.accuracy.begin(), e.accuracy.end(), e.accuracy.begin());
  Result coeff(1 / static_cast<Result>(num_examples));
  sdca_blas_scal(static_cast<blas_int>(num_classes), coeff, &e.accuracy[0]);
}


template <typename Int,
          typename Result,
          typename Context>
inline void
eval_end(
    const Int,
    const Int num_examples,
    const Context& ctx,
    eval_test<Result, multilabel_output>& e
  ) {
  // Compute the overall primal/dual objectives and their individual terms
  ctx.objective.update_primal_loss(e.primal_loss);

  e.rank_loss /= static_cast<Result>(num_examples);
}


template <typename Int,
          typename Dataset,
          typename Context>
inline void
eval_recompute_primal(
    const Int,
    const Int,
    const Dataset&,
    const Context&
  ) {
}


template <typename Int,
          typename Data,
          typename Result,
          typename Output,
          typename Context>
inline void
eval_recompute_primal(
#ifdef SDCA_ACCURATE_MATH
    const Int num_classes,
    const Int num_examples,
    const dataset<feature_input<Data>, Output, eval_train<Result, Output>>& d,
    const Context& ctx
  ) {
  // Let W = W0 + X * A'
  // (recompute W to minimize the accumulated numerical error)
  // NOTE: this should be done (if at all) for the training dataset only!
  auto D = static_cast<blas_int>(d.num_dimensions());
  auto M = static_cast<blas_int>(num_classes);
  auto N = static_cast<blas_int>(num_examples);
  if (ctx.is_prox()) {
    sdca_blas_copy(D * M, ctx.primal_initial, ctx.primal_variables);
    sdca_blas_gemm(D, M, N, d.in.features, D, ctx.dual_variables, M,
      ctx.primal_variables, CblasNoTrans, CblasTrans, 1, 1);
  } else {
    sdca_blas_gemm(D, M, N, d.in.features, D, ctx.dual_variables, M,
      ctx.primal_variables, CblasNoTrans, CblasTrans);
  }
#else
    const Int,
    const Int,
    const dataset<feature_input<Data>, Output, eval_train<Result, Output>>&,
    const Context&
  ) {
#endif
}


template <typename Int,
          typename Data,
          typename Result,
          typename Output,
          typename Context>
inline void
eval_recompute_primal(
    const Int num_classes,
    const Int num_examples,
    const dataset<model_input<Data>, Output, eval_train<Result, Output>>& d,
    const Context& ctx
  ) {
  // NOTE: this should be done always (not only with SDCA_ACCURATE_MATH),
  // but only for the training dataset.
  // (X's maintained prior to this step may be inaccurate,
  // see the update step for details)
  ctx.objective.compute_features(d.num_dimensions(), num_examples, num_classes,
    d.in.model, ctx.dual_variables, ctx.primal_initial, ctx.primal_variables);
}

}

#endif
