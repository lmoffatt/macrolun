#ifndef QM_KINETIC_MODEL_H
#define QM_KINETIC_MODEL_H

#include <qm_transition_probability.h>

template <class initial_state,class final_states>
struct  Transition_probability
{
  using type=x_i<pr<up<final_states>,dn<initial_state>>,v<double,p<>>>;
};
template <class initial_state,class final_states>
using Transition_probability_t=typename Transition_probability<initial_state,final_states>::type;

template <class initial_state,class final_states>
struct Transition_probability_rate
{
  using type=der_t<Transition_probability_t<initial_state,final_states>,v<double,p<u<s,1>>>>;
};
template <class initial_state,class final_states>
using Transition_probability_rate_t=typename Transition_probability_rate<initial_state,final_states>::type;



template <class initial_state,class final_states, class Agonist>
using Association_probability_rate_t=der_t<Transition_probability_rate_t<initial_state,final_states>,x_i<Agonist,v_t<Agonist>>>;





struct kinetic_rate{using T=double; using unit=p<u<s,-1>>;constexpr static auto className= my_static_string("kinetic_rate");};

template<class Agonist>
struct association_rate
{using T=double; using unit=decltype(p<u<s,-1>>{}/typename Agonist::unit{});constexpr static auto className=
      Agonist::className+my_static_string("association_rate");};


template<class A, std::size_t N>
struct r
{
  constexpr auto n(A){ return N;}
  constexpr static auto className=A::className+to_static_string<N>();
  using Ag=A;
  constexpr auto diff(std::size_t N2) {if (N>N2) return N-N2; else return 0ul;}

};






template <std::size_t N, class...As>
struct S: As...
{ using As::n...;
  constexpr static auto className=my_static_string("S")+(to_static_string<N>()+...+As::className);
  constexpr auto operator[](std::integral_constant<std::size_t,N>)const {return *this;}

};




struct Valid{};
struct Invalid{};





template <std::size_t N1, class...A1s, std::size_t N2,  class...A2s>
constexpr auto diff(S<N1,A1s...> one,S<N2,A2s...> )
{
  return (Valid{}|...|[&one](auto v){
    if constexpr (std::is_same_v<v, Invalid>) return v;
    else
        {
      auto n=A2s{}.diff(one.n(typename A2s::Ag{}));
      if constexpr (n==0) return v;
      else if constexpr (n>1) return Invalid{};
      else if constexpr (std::is_same_v<v,Valid>) return typename A2s::Ag{};
      else return Invalid{};
    }
  });
}


template<class...> struct k_s;




template <std::size_t N1, class...A1s, std::size_t N2>
struct k_s<S<N1,A1s...>, S<N2,A1s...>> {using type=Transition_probability_rate_t<S<N1,A1s...>,S<N2,A1s...>>;};

template <std::size_t N1, class...A1s, std::size_t N2, class...A2s>
struct k_s<S<N1,A1s...>, S<N2,A2s...>> {

  using res=decltype (kind(S<N1,A1s...>{},S<N2,A2s...>{}));

  using type=std::conditional_t<std::is_same_v<Invalid,res>,Nothing,
      std::conditional_t<std::is_same_v<Valid, res>,Transition_probability_rate_t<S<N1,A1s...>,S<N2,A1s...>>,
                                                     Association_probability_rate_t<S<N1,A1s...>,S<N2,A2s...>,res>>>;
};






template<std::size_t I, std::size_t...Is>
struct Li
{
  constexpr auto operator[](std::integral_constant<std::size_t,I>)const {return *this;}

};

template<class ...Ls>
struct Con:Ls...
{
  using Ls::operator[]...;
  template<std::size_t I>
  constexpr auto operator()()const
  {
    return (*this)[std::integral_constant<std::size_t,I>{}];
  }
};


template <class... Ss>
struct States:Ss...
{
  using Ss::operator[]...;
  template<std::size_t I>
  constexpr auto operator()()const
  {
    return (*this)[std::integral_constant<std::size_t,I>{}];

  }
};


template<class, class > struct StateModel;



template <class... Ss,class ...Ls>
struct StateModel<States<Ss...>,Con<Ls...>>
{
  template<std::size_t I>
  using S_t=decltype (States<Ss...>{}[std::integral_constant<std::size_t,I>{}]);

  template<std::size_t I, std::size_t J>
  using k_t=typename k_s<S_t<I>,S_t<J>>::type;




};






template <class initial_state,class final_states>
using k_t=Transition_probability_rate_t<initial_state,final_states>;


template <class initial_state,class final_states, class Agonist>
using ka_t=Association_probability_rate_t<initial_state,final_states, Agonist>;






struct Q0{constexpr static auto className=my_static_string("Q0");};
struct Qa{constexpr static auto className=my_static_string("Qa");};





typedef Cs<> param;








#endif // QM_KINETIC_MODEL_H
