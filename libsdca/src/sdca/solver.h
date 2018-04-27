#ifndef SDCA_SOLVER_H
#define SDCA_SOLVER_H

#include <cassert>
#include <numeric>
#include <random>

#include "solver/context.h"
#include "solver/eval.h"
#include "solver/update.h"

namespace sdca {

template <typename Data,
          typename Result,
          template <typename> class Input,
          typename Output,
          template <typename, typename> class Objective>
class solver {
public:
  typedef Data data_type;
  typedef Result result_type;
  typedef Input<Data> input_type;
  typedef Output output_type;
  typedef Objective<Data, Result> objective_type;

  typedef solver_context<Data, Result, Input, Output, Objective> context_type;
  typedef solver_scratch<Data, Input> scratch_type;


  explicit solver(
      context_type& __context
    ) :
      ctx_(__context),
      is_evaluated_(false)
  {}


  void solve() {
    begin_solve();
    while (ctx_.status == solver_status::solving) {

      begin_epoch();
      for (auto& example : examples_) {
        update_variables(example, ctx_, scratch_);
      }

      end_epoch();
    }
    end_solve();
  }


protected:
  context_type& ctx_;
  scratch_type scratch_;

  bool is_evaluated_;
  std::minstd_rand generator_;
  std::vector<size_type> examples_;


  void begin_solve() {
    reporting::begin_solve(ctx_);

    ctx_.status = (ctx_.criteria.max_epoch > ctx_.epoch)
                  ? solver_status::solving
                  : solver_status::max_epoch;

    scratch_.init(ctx_.train);
    if (ctx_.criteria.eval_on_start) {
      evaluate_solution();
      check_stopping_criteria<Data, Result>(ctx_);
    }

    if (ctx_.status == solver_status::solving) {
      generator_.seed();
      examples_.resize(ctx_.train.num_examples());
      std::iota(examples_.begin(), examples_.end(), 0);

      ctx_.solve_time.resume();
    }
  }


  void end_solve() {
    ctx_.solve_time.stop();
    if (!is_evaluated_) {
      evaluate_solution();
    }

    reporting::end_solve(ctx_);
  }


  void begin_epoch() {
    is_evaluated_ = false;
    if (need_shuffle<input_type>::value)
      std::shuffle(examples_.begin(), examples_.end(), generator_);
  }


  void end_epoch() {
    ctx_.solve_time.stop();

    ++ctx_.epoch;
    if ((ctx_.criteria.eval_epoch > 0) &&
        (ctx_.epoch % ctx_.criteria.eval_epoch == 0)) {
      evaluate_solution();
    }

    check_stopping_criteria<Data, Result>(ctx_);

    reporting::end_epoch(ctx_, is_evaluated_);

    ctx_.solve_time.resume();
  }


  void evaluate_solution() {
    ctx_.eval_time.resume();

    evaluate_dataset(ctx_, ctx_.train, scratch_);

    size_type id(0);
    for (auto& test_set : ctx_.test) {
      evaluate_dataset(ctx_, test_set, scratch_, id++);
    }

    is_evaluated_ = true;
    ctx_.eval_time.stop();
  }

};


template <typename Data,
          typename Result,
          template <typename> class Input,
          typename Output,
          template <typename, typename> class Objective>
inline solver<Data, Result, Input, Output, Objective>
make_solver(
    solver_context<Data, Result, Input, Output, Objective>& ctx
  ) {
  return solver<Data, Result, Input, Output, Objective>(ctx);
}

}

#endif
