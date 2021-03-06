#include "FeatureCache.hh"

#include <string.h>

using namespace Ami;

FeatureCache::FeatureCache() :
  _update(false)
{
}

FeatureCache::~FeatureCache()
{
}

void     FeatureCache::clear()
{
  _names.clear();
  _cache.clear();
  _damaged.clear();
}

unsigned FeatureCache::add(const std::string& name)
{ return add(name.c_str()); }

unsigned FeatureCache::add(const char* name)
{
  std::string sname(name);
  
  for(unsigned k=0; k<_names.size(); k++)
    if (sname==_names[k])
      return k;

  _update=true;
  _names.push_back(sname);
  unsigned len = _names.size();
  unsigned index = len - 1;
  _cache  .resize(len);
  _damaged.resize((len+0x1f)>>5);
  _damaged[index>>5] |= 1<<(index&0x1f);
  return index;
}

void FeatureCache::add(const FeatureCache& c)
{
  for(unsigned k=0; k<c._names.size(); k++)
    add(c._names[k]);
}

int         FeatureCache::lookup(const char* name) const
{
  std::string sname(name);

  for(unsigned k=0; k<_names.size(); k++)
    if (sname==_names[k])
      return int(k);

  return -1;
}

unsigned    FeatureCache::entries() const { return _names.size(); }
const std::vector<std::string>& FeatureCache::names  () const { return _names; }
double      FeatureCache::cache  (int index, bool* damaged) const 
{
  if (damaged) *damaged = (index<0) || ((_damaged[index>>5]>>(index&0x1f)) & 1);
  return index>=0 ? _cache[index] : 0;
}

void        FeatureCache::cache  (int index, double v, bool damaged)
{
  if (index>=0) {
    _cache[index] = v;
    uint32_t mask = 1<<(index&0x1f);
    if (damaged)
      _damaged[index>>5] |= mask;
    else
      _damaged[index>>5] &= ~mask;
  }
}

void        FeatureCache::cache  (const FeatureCache& c)
{
  for(unsigned i=0,j=0; i<c._cache.size(); i++,j++) {
    while(_names[j]!=c._names[i]) j++;
    cache(j, c._cache[i], (c._damaged[i>>5] & (1<<(i&0x1f))));
  }
}

char*  FeatureCache::serialize(int& len) const
{
  len = _names.size()*FEATURE_NAMELEN;
  char* result = new char[len];
  char* p = result;
  for(unsigned k=0; k<_names.size(); k++, p+=FEATURE_NAMELEN) {
    strncpy(p, _names[k].c_str(), FEATURE_NAMELEN);
  }
  return result;
}

bool   FeatureCache::update()
{
  bool upd(_update);
  _update=false;
  return upd;
}

void   FeatureCache::dump() const
{
  printf("FeatureCache entries %d\n",(int) _names.size());
  for(unsigned k=0; k<_names.size(); k++) {
    printf("  %s\n",_names[k].c_str());
  }
}

void   FeatureCache::start()
{
  for(unsigned k=0; k<_damaged.size(); k++)
    _damaged[k] = -1UL;
}
