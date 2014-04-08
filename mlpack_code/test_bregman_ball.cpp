#include "bregman_ball.hpp"

#include "KLDivergence.hpp"
#include "L2Divergence.hpp"

using namespace bmst;

int main (int argc, char* argv[])
{
  std::vector<double> q_vec(5, 0.1);
  q_vec[0] = 0.5;
  Point<double> q(q_vec);

  std::vector<double> mu_vec(5,0.2);  
  Point<double> mu(mu_vec);
  
  //d (q, mu) = 0.18088649371309945
  
  std::cout << "Testing KL ball.\n";
  
  BregmanBall<double, KLDivergence<double> > kl_ball_small(mu, 0.001);
  BregmanBall<double, KLDivergence<double> > kl_ball_large(mu, 0.1);

  // y has a neighbor within 0.001
  bool can_prune_small = kl_ball_small.CanPruneRight(q, 0.05);
  //std::cout << "small prune: " << can_prune_small << "\n";
  assert(can_prune_small);
  
  bool can_prune_large = kl_ball_large.CanPruneRight(q, 0.05);
  //std::cout << "large prune: " << can_prune_large << "\n";
  assert(!can_prune_large);
  
  std::cout << "KL Ball Passed.\n";
  
  
  std::cout << "Testing L2 Ball\n";
  
  // 0.065
  
  BregmanBall<double, L2Divergence<double> > l2_ball_small(mu, 0.001);
  BregmanBall<double, L2Divergence<double> > l2_ball_large(mu, 10.0);
  
  can_prune_small = l2_ball_small.CanPruneRight(q, 0.03);
  //std::cout << "small prune: " << can_prune_small << "\n";
  assert(can_prune_small);
  
  can_prune_large = l2_ball_large.CanPruneRight(q, 0.03);
  //std::cout << "large prune: " << can_prune_large << "\n";
  assert(!can_prune_large);
  
  std::cout << "L2 Ball Passed.\n";
  
  return 0;
}
