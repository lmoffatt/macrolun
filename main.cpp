#include <iostream>
#include <qm_tensor_model.h>
#include <qm_transition_probability.h>
#include <qm_kinetic_model.h>

struct ms{constexpr static auto  className=my_static_string("ms");};

struct ns{constexpr static auto  className=my_static_string("nsample");};
struct uMATP{constexpr static auto  className=my_static_string("uM_ATP");};
struct pA{constexpr static auto className=my_static_string("pA");};

using ms_u=p<u<ms,1>>;
using ns_u=p<u<ns,1>>;
using uMATP_u=p<u<uMATP,1>>;
using pA_u=p<u<pA,1>>;


struct mytime
{
  typedef double T;
  typedef ms_u unit;
  constexpr static auto  className=my_static_string("time");
};
struct num_samples
{
  typedef std::size_t T;
  using unit=ns_u;
  constexpr static auto  className=my_static_string("num_samples");
};

struct con_ATP {using T=double; using unit=uMATP_u; constexpr static auto className=my_static_string("con_ATP"); };

struct i_Patch {using T=double; using unit=pA_u; constexpr static auto className= my_static_string("i_Patch");};



int main(int argc, char **argv)
{
  std::cerr<<argv[0]<<"\n";
  std::cerr<<argv[1]<<"\n";


  auto priors=D(param{},Normal_Distribution{},mean<kinetic_rate>{},stddev<kinetic_rate>{});




  auto q=quimulun
  {
          x_i(mytime{},vec<mytime>{}),
          x_i(num_samples{},vec<mytime>{}),
          x_i(con_ATP{},vec<mytime>{})




      };




};
