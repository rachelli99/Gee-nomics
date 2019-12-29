#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Trie.h"
#include <algorithm>
#include <unordered_map>
using namespace std;

bool waysToCompare(GenomeMatch x, GenomeMatch y) {return (x.percentMatch > y.percentMatch);}

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength;
    Trie<pair<int, int>> m_trie;
    int genomeIndex; //total genome number
    vector<Genome> vg;
    
    vector<DNAMatch> setMapToVector(unordered_map<string, pair<int, int>> m1, vector<DNAMatch> v1) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSearchLength = minSearchLength;
    genomeIndex = 0;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    vg.push_back(genome);
    if(genome.length() < m_minSearchLength) return;
    
    for(int x = 0; x <= genome.length() - minimumSearchLength(); x++)
    {
        string extracted  = "";
        genome.extract(x, m_minSearchLength, extracted);
        pair<int, int> value;
        value.first = genomeIndex;
        value.second = x;
        m_trie.insert(extracted, value);
    }
    genomeIndex++;
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    if(fragment.length() < minimumLength || minimumLength < minimumSearchLength()) return false;
    vector<pair<int, int>> result;
    unordered_map<string, pair<int, int>> map;
    
    //check the first fragment from 0 to minimumSearchLength to see if there's match
    string search = fragment.substr(0, minimumSearchLength());

    result = m_trie.find(search, exactMatchOnly);

    if(result.empty()) return false;
    
    vector<pair<int, int>>::iterator it = result.begin();
    for(; it != result.end(); it++)
    {
        int index = (*it).first;
        int pos = (*it).second;
        Genome g = vg[index];  //retrive the Genome from the vector
        
        string segment;
        g.extract(pos, (int)fragment.length(), segment);
        
        int misMatches = 0;
        int x;
        for(x = 0; x < fragment.length(); x++)
        {
            if(exactMatchOnly)
            {
                if(segment[x] != fragment[x])
                    break;
            }
            else
            {
                if(segment[x] != fragment[x]) misMatches++;
                if(misMatches == 2) break; //single mismatch allowed
            }
        }
        
        int matchLength = (int)segment.substr(0, x).length();
        if(matchLength >= minimumLength)
        {
            if(map.find(g.name()) == map.end()) //genome fragment doesn't exist in map
            {
                //insert the current genome into the map in the format of <pos, length>
                pair<int, int> tempPair (pos, matchLength);
                pair<string, pair<int, int>> tempMap (g.name(), tempPair);
                map.insert(tempMap);
            }
            else //genome exist, compare length
            {
                if(map[g.name()].second < matchLength)
                {
                    map[g.name()].first = pos;
                    map[g.name()].second = matchLength;
                }
            }
        }
    }
    matches = setMapToVector(map, matches);
    
    if(matches.empty()) return false;
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if(fragmentMatchLength < minimumSearchLength()) return false;
    unordered_map<string, int> map; //mapping name to count
    
    for(int x = 0; x < query.length(); x += fragmentMatchLength)
    {
        //extract matching fragment and push into dna vector
        string fragment;
        query.extract(x, fragmentMatchLength, fragment);
        vector<DNAMatch> dna;
        findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, dna);
        
        //set non-duplicated DNA and their occurance into the map
        vector<DNAMatch>::iterator it = dna.begin();
        for(; it != dna.end(); it++)
        {
            if(map.find((*it).genomeName) == map.end()) //no such genome
            {
                pair<string, int> tempPair ((*it).genomeName, 1);
                map.insert(tempPair);
            }
            else //genome exist
            {
                map[(*it).genomeName]++; //increament the second integer count value
            }
        }
    }
    
    int segment = query.length() / fragmentMatchLength;
    //iterate through the map to find the amount above threshold to results
    unordered_map<string, int>::iterator it2 = map.begin();
    for (; it2 != map.end(); it2++)
    {
        int numMatch = (*it2).second;
        double percentage = 100*numMatch/segment;
        if(percentage >= matchPercentThreshold)
        {
            GenomeMatch gm;
            gm.genomeName = (*it2).first;
            gm.percentMatch = percentage;
            results.push_back(gm);
        }
    }
    
    //sorts in decending order
    sort(results.begin(), results.end(), waysToCompare);
    
    if(results.empty()) return false;
    return true;
}

//-----------------------------------private function--------------------------------------

vector<DNAMatch> GenomeMatcherImpl::setMapToVector(unordered_map<string, pair<int, int>> m1, vector<DNAMatch> v1) const
{
    unordered_map<string, pair<int, int>>::iterator it = m1.begin();
    for(; it != m1.end(); it++)
    {
        DNAMatch d;
        d.genomeName = (*it).first;
        d.position = (*it).second.first;
        d.length = (*it).second.second;
        v1.push_back(d);
    }
    return v1;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
