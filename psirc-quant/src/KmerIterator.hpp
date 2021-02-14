#ifndef BFG_KMER_ITERATOR_HPP
#define BFG_KMER_ITERATOR_HPP

#include <iterator>
#include "Kmer.hpp"


/* Short description:
 *  - Easily iterate through kmers in a read
 *  - If the read contains any N, then the N is skipped and checked whether
 *    there is a kmer to the right of the N
 * */
class KmerIterator : public std::iterator<std::input_iterator_tag, std::pair<Kmer, int>, int> {
 public:
  KmerIterator() : s_(NULL), p_(), invalid_(true){}
  KmerIterator(const char *s) : s_(s), p_(), invalid_(false){ find_next(-1,-1,false);}
  KmerIterator(const KmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_) {}

  KmerIterator& operator++();
  KmerIterator operator++(int);
  void raise(Kmer& km, Kmer& rep);
  void jumpTo(int pos);

  bool operator==(const KmerIterator& o);
  bool operator!=(const KmerIterator& o) { return !this->operator==(o);}

  std::pair<Kmer, int>& operator*();
  std::pair<Kmer, int> *operator->();

 protected:
  void find_next(size_t i, size_t j, bool last_valid);

  const char *s_;
  std::pair<Kmer, int> p_;//p_->first is the k-mer in bits representation; p_->second is the start position of the kmer
  bool invalid_;
  //bool isLinear_;//whether the sequence is linear(true) or circular(false);
};

/*
 * Extension from KmerIterator for iterating through kmers in both linear and circular trancript(form a circle in de Bruijn graph). 
 * */
class TransKmerIterator : public KmerIterator{
 public:
  TransKmerIterator(){}
  TransKmerIterator(const char* s, size_t len, bool isL=true) : KmerIterator(s), len_(len), isLinear_(isL){}
  //TransKmerIterator(const TransKmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_), len_(o.len_), isLinear_(o.isLinear_) {}
  TransKmerIterator(const TransKmerIterator& o) : KmerIterator(o), len_(o.len_), isLinear_(o.isLinear_) {}
  
  TransKmerIterator& operator++();
  TransKmerIterator operator++(int);
  void raise(Kmer& km, Kmer& rep);
  void jumpTo(int pos);
  
  bool operator==(const TransKmerIterator& o);
  bool operator!=(const TransKmerIterator& o) { return !this->operator==(o);}
 protected:
  void find_next(size_t i, size_t j, bool last_valid);//find next valid k-mers from kmer in [i, j]

  size_t len_;
  bool isLinear_;
};
#endif // BFG_KMER_ITERATOR_HPP
