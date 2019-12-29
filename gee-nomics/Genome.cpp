#include "provided.h" 
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
    bool operator<(const Genome& other) const {return (m_name < other.name());}
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    vector<Genome> result;
    
    if(!genomeSource) return false;
    
    string line;
    string name = "";
    string genes = "";
    while(getline(genomeSource, line))
    {
        if(name == "" && genes != "") return false; //no name but there's genes
        
        if(line[0] == '>')
        {
            if(line.length() > 1)
            {
                if(name != "") //already filled with a complete genome
                {
                    Genome g(name, genes);
                    result.push_back(g);
                    name = "";
                    genes = "";
                }
                name = line.substr(1);
                continue;
            }
            else return false; //empty name line
        }
        else //genes line
        {
            for(int x = 0; x < line.length(); x++)
            {
                toupper(line[x]); //make all uppercase
            }
            for(int x = 0; x < line.length(); x++) //all non name character are A,C,T,G,N
            {
                if(line[x] != 'A' && line[x] != 'C' && line[x] != 'T' && line[x] != 'G' && line[x] != 'N')
                {
                    return false;
                }
            }
            genes += line;
        }
    }
    
    if(name != "" && genes != "") //push back the last genome stored in name and genes
    {
        Genome g(name, genes);
        result.push_back(g);
    }
    
    genomes = result;
    return true;
}

int GenomeImpl::length() const
{
    return (int)m_sequence.length();
}

string GenomeImpl::name() const
{
    return m_name;  // This compiles, but may not be correct
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if((position + length) > m_sequence.length()) return false;
    
    fragment = m_sequence.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
