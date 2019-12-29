#include "Trie.h"
#include <iostream>
#include <vector>
#include "provided.h"
#include <fstream>
#include <cassert>
using namespace std;

void test_loadAndExtract();
void test_trie();
void test_extract();
void test_GenomeImpl_Exceptions();
void test_GenomeMatcher1();
void test_GenomeMatcher2();
void test_GenomeMatcher3();
void test_GenomeMatcher();

//int main(int argc, const char * argv[])
//{
//    test_GenomeMatcher();
//}

void test_loadAndExtract()
{
    //string filename = "c:/genomes/Ferroplasma_acidarmanus.txt";
    string filename = "/Users/Rachel/Desktop/data/Halorientalis_regularis.txt";
    // Open the data file and get a ifstream object that can be used to read its
    // contents.
    ifstream strm(filename);
    if (!strm)
    {
        cout << "Cannot open " << filename << endl;
    }
    vector<Genome> vg;
    bool success = Genome::load(strm, vg); // Load the data via the stream.
    if (success)
    {
        cout << "Loaded " << vg.size() << " genomes successfully:" << endl;
        for (int k = 0; k != vg.size(); k++)
        {
            cout << vg[k].name() << endl;
            cout << vg[k].length() << endl;
        }
    }
    else
        cout << "Error loading genome data" << endl;
    
    Genome g = vg[0];
    string fragment = "";
    g.extract(1000, 6, fragment);
    cout << fragment << endl;
    
}

void test_trie()
{
    Trie<int> trie; // This is like std::multimap<std::string,int> trie;
    //    trie.insert("hit", 1);
    //    trie.insert("hit", 2);
    //    trie.insert("hip", 4);
    trie.insert("GATTACA",42); //GATTACAà{42} //
    trie.insert("GATTACA",17); //GATTACAà{42,17}
    trie.insert("GATTACA", 42); // GATTACA à {42, 17, 42}
    trie.insert("GCTTACA",30); //GCTTACAà{30}
    //    trie.print();
    vector<int> v = trie.find("GATTACA", false);
    for(int x = 0; x != v.size(); x++)
        cout << v[x] << "  " << endl;
}

void test_extract()
{
    Genome g("oryx", "GCTCGGNACACATCCGCCGCGGACGGGACGGGATTCGGGCTGTCGATTGTCTCACAGATCGTCGACGTACATGACTGGGA");
    string f1, f2, f3;
    bool result1 = g.extract(0, 5, f1); // result1 = true, f1 = “GCTCG”;
    bool result2 = g.extract(74, 6, f2); // result2 = true, f2 = “CTGGGA”;
    bool result3 = g.extract(74, 7, f3); // result3 = false, f3 is unchanged
    
    assert(g.length()==80 && g.name()=="oryx");
    assert(result1 && f1=="GCTCG");
    assert(result2 && f2=="CTGGGA");
    assert(!result3 && f3=="");
}

void test_GenomeImpl_Exceptions()
{
    vector<Genome> v;
    ifstream inputFile("/Users/Rachel/Desktop/data/Halorientalis_regularis.txt");
    assert(!Genome::load(inputFile, v));
    
    ifstream inputFile2("/Users/Rachel/Desktop/data/Halorientalis_regularis.txt");
    assert(!Genome::load(inputFile2, v));
    
    ifstream inputFile3("/Users/Rachel/Desktop/data/Halorientalis_regularis.txt");
    assert(!Genome::load(inputFile3, v));
    
    ifstream inputFile4("/Users/Rachel/Desktop/data/Halorientalis_regularis.txt");
    assert(!Genome::load(inputFile4, v));
    assert(v.empty());
    
    ifstream inputFile5("/Users/Rachel/Desktop/data/Halorientalis_regularis.txt");
    assert(!Genome::load(inputFile5, v));
    
}

void test_GenomeMatcher1()
{
    Genome g1("Genome_1","CGGTGTACNACGACTGGGGATAGAATATCTTGACGTCGTACCGGTTGTAGTCGTTCGACCGAAGGGTTCCGCGCCAGTAC");
    Genome g2("Genome_2","TAACAGAGCGGTNATATTGTTACGAATCACGTGCGAGACTTAGAGCCAGAATATGAAGTAGTGATTCAGCAACCAAGCGG");
    Genome g3("Genome_3","TTTTGAGCCAGCGACGCGGCTTGCTTAACGAAGCGGAAGAGTAGGTTGGACACATTNGGCGGCACAGCGCTTTTGAGCCA");
    
    
    GenomeMatcher gm1(4);
    gm1.addGenome(g1);
    gm1.addGenome(g2);
    gm1.addGenome(g3);
    
    std::vector<DNAMatch> matches;
    bool result;
    result = gm1.findGenomesWithThisDNA("GAAG", 4, true, matches);
    assert(result);
    int SIZE=matches.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches[i].genomeName<<" of length "<<matches[i].length<<" at position "<<matches[i].position<<endl;
    }
    cout<<endl;
    
    
    std::vector<DNAMatch> matches_2;
    result = gm1.findGenomesWithThisDNA("GAATAC", 4, true, matches_2);
    SIZE=matches_2.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_2[i].genomeName<<" of length "<<matches_2[i].length<<" at position "<<matches_2[i].position<<endl;
    }
    assert(result);
    cout << endl;
    
    std::vector<DNAMatch> matches_3;
    result = gm1.findGenomesWithThisDNA("GAATAC", 6, true, matches_3);
    assert(!result);
    cout<<endl;
    
    
    
    std::vector<DNAMatch> matches_4;
    result = gm1.findGenomesWithThisDNA("GAATAC", 6, false, matches_4);
    SIZE=matches_4.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_4[i].genomeName<<" of length "<<matches_4[i].length<<" at position "<<matches_4[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    cout<<"Test 5: GTATAT, 6, false, matches:"<<endl;
    std::vector<DNAMatch> matches_5;
    result = gm1.findGenomesWithThisDNA("GTATAT", 6, false, matches_5);
    SIZE=matches_5.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_5[i].genomeName<<" of length "<<matches_5[i].length<<" at position "<<matches_5[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    cout<<"Test 6: GAATACG, 6, false, matches:"<<endl;
    std::vector<DNAMatch> matches_6;
    result = gm1.findGenomesWithThisDNA("GAATACG", 6, false, matches_6);
    SIZE=matches_6.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_6[i].genomeName<<" of length "<<matches_6[i].length<<" at position "<<matches_6[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    result=false;
    cout<<"Test 7: GAAGGGTT, 5, false, matches:"<<endl;
    std::vector<DNAMatch> matches_7;
    result = gm1.findGenomesWithThisDNA("GAAGGGTT", 5, false, matches_7);
    SIZE=matches_7.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_7[i].genomeName<<" of length "<<matches_7[i].length<<" at position "<<matches_7[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    result=false;
    cout<<"Test 8: GAAGGGTT, 6, false, matches:"<<endl;
    std::vector<DNAMatch> matches_8;
    result = gm1.findGenomesWithThisDNA("GAAGGGTT", 6, false, matches_8);
    SIZE=matches_8.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_8[i].genomeName<<" of length "<<matches_8[i].length<<" at position "<<matches_8[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    result=false;
    cout<<"Test 9: ACGTGCGAGACTTAGAGCC, 12, false, matches:"<<endl;
    std::vector<DNAMatch> matches_9;
    result = gm1.findGenomesWithThisDNA("ACGTGCGAGACTTAGAGCC", 12, false, matches_9);
    SIZE=matches_9.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_9[i].genomeName<<" of length "<<matches_9[i].length<<" at position "<<matches_9[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    
    result=false;
    cout<<"Test 10: ACGTGCGAGACTTAGAGCG, 12, false, matches:"<<endl;
    std::vector<DNAMatch> matches_10;
    result = gm1.findGenomesWithThisDNA("ACGTGCGAGACTTAGAGCG", 12, false, matches_10);
    SIZE=matches_10.size();
    for(int i=0;i<SIZE;i++)
    {
        cout<<matches_10[i].genomeName<<" of length "<<matches_10[i].length<<" at position "<<matches_10[i].position<<endl;
    }
    cout<<endl;
    assert(result);
    
    std::vector<DNAMatch> matches_11;
    result = gm1.findGenomesWithThisDNA("GAAG", 3, true, matches_11);
    assert(!result);
    result=true;
    result = gm1.findGenomesWithThisDNA("GAAG", 5, true, matches_11);
    assert(!result);
}

void test_GenomeMatcher2()
{
    Genome g1("Genome_a","ACGTATACNACGAGCGGGGATAGAATTTTCTGACGTCGTACCGGTTGTAGTCGTTCGACCGAAGGGTTCCGCGCCAGTAC");
    Genome g2("Genome_b","TAACAGAGCGGTNATATTGTTACGAATCACGTGCGAGACTTAGAGCCAGAATATGAAGTAGTGATTCAGCAACCAAGCGG");
    Genome g3("Genome_c","TTTTGAGCCAGCGACGCGGCTTGCTTAACGAAGCGGAAGAGTAGGTTGGACACATTNGGCGGCACAGCGCTTTTGAGCCA");
    //                    *1   *6   *11  *16  *21  *26  *31  *36  *41  *46  *51  *56  *61  *66  *71  *76
    
    Genome g4("Genome_A","ACGTACNTNCAA");
    Genome g5("Genome_B","ACGAACGTACAAANAA");
    Genome g6("Genome_C","ACGTACGTACGT");
    
    Genome g7("Genome_C","ACGTACGTACGT");
    Genome g8("Genome_C","ACGT");
    Genome g9("Genome_C","ACGTACGTACGT");
    
    GenomeMatcher gm1(4);
    gm1.addGenome(g1);
    gm1.addGenome(g2);
    gm1.addGenome(g3);
    gm1.addGenome(g4);
    gm1.addGenome(g5);
    gm1.addGenome(g6);
    
    std::vector<GenomeMatch> results;
    Genome query("Query","ACGTACNTANTTTTC");
    gm1.findRelatedGenomes(query, 5, true, 30, results);
    
    int SIZE=results.size();
    cout.setf(ios::fixed);
    cout.precision(2);
    for(int i=0;i<SIZE;i++)
    {
        cout<<results[i].genomeName<<" has percentage: "<<results[i].percentMatch<<endl;
    }
}

void test_GenomeMatcher3()
{
    Genome g5("Genome_sameas_a_2","AAAAGTCCCCTTTTATCCCA");
    Genome g1("Genome_a","AAAAGTCCCCTTTTATCCCA");
    Genome g2("Genome_b","TCCCNGAAATTTTTGCCCCCNNNNN");
    Genome g3("Genome_c","TTTTGAGCCAGCGACGCGGCTTGCTTAACGAAGCGGAAGAGTAGGTTGGACACATTNGGCGGCACAGCGCTTTTGAGCCA");
    Genome g4("Genome_sameas_a","AAAAGTCCCCTTTTATCCCA");
    
    //                    *1   *6   *11  *16  *21  *26  *31  *36  *41  *46  *51  *56  *61  *66  *71  *76
    GenomeMatcher gm1(4);
    gm1.addGenome(g1);
    gm1.addGenome(g2);
    gm1.addGenome(g3);
    gm1.addGenome(g4);
    gm1.addGenome(g5);
    //    gm1.addGenome(g6);
    
    std::vector<GenomeMatch> results;
    Genome query("Query","AAAAGCNNNATTTTGACGNTCCCNCATATATCCCAGAAATCCCCC");
    //                    *1   *6   *11  *16  *21  *26  *31  *36  *41
    
    gm1.findRelatedGenomes(query, 4, true, 11.11, results);
    int SIZE=results.size();
    cerr << "result size: " << results.size() << endl;
    cout.setf(ios::fixed);
    cout.precision(2);
    for(int i=0;i<SIZE;i++)
    {
        cout<<results[i].genomeName<<" has percentage: "<<results[i].percentMatch<<endl;
    }
}

void test_GenomeMatcher()
{
    Genome g3("Genome 3", "TTTTGAGCCAGCGACGCGGCTTGCTTAACGAAGCGGAAGAGTAGGTTGGACACATTNGGCGGCACAGCGCTTTTGAGCCA");
    Genome g2("Genome 2", "TAACAGAGCGGTNATATTGTTACGAATCACGTGCGAGACTTAGAGCCAGAATATGAAGTAGTGATTCAGCAACCAAGCGG");
    Genome g1("Genome 1", "CGGTGTACNACGACTGGGGATAGAATATCTTGACGTCGTACCGGTTGTAGTCGTTCGACCGAAGGGTTCCGCGCCAGTAC");
    
//    Genome query("query", "ACTGACTGTCGATTGACATGCATGGCTA");
    Genome query("query", "CGGTGTACNACGACTGGGGATAGAATATCTTGACGTCGTAC");
    
    GenomeMatcher gm1(4);
    gm1.addGenome(g1);
    gm1.addGenome(g2);
    gm1.addGenome(g3);
    std::vector<DNAMatch> matches;
    std::vector<GenomeMatch> results;
    bool result;
//    result = gm1.findGenomesWithThisDNA("GAAG", 5, true, matches);
//    cout << result <<endl;
//    for(int x = 0; x < matches.size(); x++)
//    {
//        cout << matches[x].genomeName << " of length " << matches[x].length<< " at position " << matches[x].position << endl;
//    }
    
    result = gm1.findRelatedGenomes(query, 5, false, 2, results);
    cout << result << endl;
}
