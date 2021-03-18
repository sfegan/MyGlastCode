//-*-mode:c++; mode:font-lock;-*-

#ifndef BLOCKPARTITION_HPP
#define BLOCKPARTITION_HPP

#include <vector>
#include <algorithm>

template<typename T> class BlockCostFunction
{
public:
  virtual ~BlockCostFunction();
  virtual double cost(typename std::vector<T>::const_iterator begin,
		      typename std::vector<T>::const_iterator end) = 0;
};

template<typename T> BlockCostFunction<T>::~BlockCostFunction()
{
  // nothing to see here
}

template<typename T> class BlockPenalty: public BlockCostFunction<T>
{
public:
  BlockPenalty(BlockCostFunction<T>& base_cost_fn, double penalty = 0):
    BlockCostFunction<T>(), m_base_cost_fn(base_cost_fn), m_penalty(penalty) {}
  virtual ~BlockPenalty();
  virtual double cost(typename std::vector<T>::const_iterator begin,
		      typename std::vector<T>::const_iterator end);
  void setPenalty(double penalty) { m_penalty=penalty; }
  double penalty() const { return m_penalty; }
private:
  BlockCostFunction<T>& m_base_cost_fn;
  double                m_penalty;
};

template<typename T> BlockPenalty<T>::~BlockPenalty()
{
  // nothing to see here
}

template<typename T> double BlockPenalty<T>::
cost(typename std::vector<T>::const_iterator begin,
     typename std::vector<T>::const_iterator end)
{
  return m_base_cost_fn.cost(begin,end) + m_penalty;
}

template<typename T> class BlockPartionCalculator
{
public:
  BlockPartionCalculator(const std::vector<T>& events,
			 BlockCostFunction<T>& cost_fn):
    m_events(events), m_cost_fn(cost_fn) { }
  ~BlockPartionCalculator();
  double compute(std::vector<unsigned>& partition);

private:
  const std::vector<T>& m_events;
  BlockCostFunction<T>& m_cost_fn;
};

template<typename T> BlockPartionCalculator<T>::~BlockPartionCalculator()
{
  // nothing to see here
}

template<typename T>
double BlockPartionCalculator<T>::compute(std::vector<unsigned>& partition)
{
  const unsigned nevents = m_events.size();
  std::vector<double> opt;
  opt.reserve(nevents+1);

  std::vector<unsigned> lastchange;
  lastchange.reserve(nevents+1);
  
  opt.push_back(m_cost_fn.cost(m_events.begin(), m_events.begin()+1));
  lastchange.push_back(0);

  for(unsigned in=1;in<nevents;in++)
    {
      double copt = m_cost_fn.cost(m_events.begin(), m_events.begin()+in+1);
      double jopt = 0;
      //std::cout << in << ' ' << 0 << ' ' << copt << '\n';
      for(unsigned ij=1;ij<=in;ij++)
	{
	  double x = 
	    m_cost_fn.cost(m_events.begin()+ij,m_events.begin()+in+1)
	    + opt[ij-1];
	  if(x<copt)copt=x,jopt=ij;
	  //std::cout << in << ' ' << ij << ' ' << x << ' ' << jopt << ' ' << copt << '\n';
	}
      opt.push_back(copt);
      lastchange.push_back(jopt);
      //std::cout << '\n';
    }
  
  partition.clear();
  for(unsigned ip = lastchange.back(); ip!=0; ip=lastchange[ip-1])
    partition.push_back(ip);
  std::reverse(partition.begin(),partition.end());

  return opt.back();
}

#endif // defined BLOCKPARTITION
