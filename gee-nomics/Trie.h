#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <iostream> //not here

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    
    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
    struct Node
    {
        std::vector<ValueType> val;
        char label;
        std::vector<Node*> children;
    };
    
    Node* root;
    
    void clearUp(Node* cur);
    std::vector<ValueType> combineVectors(std::vector<ValueType> v1, std::vector<ValueType> v2) const;
};

//---------------Trie.cpp---------------
template<typename ValueType>
Trie<ValueType>::Trie()
{
    root = new Node;
    
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
    clearUp(root);
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    clearUp(root);
    root = new Node;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
    if(key.size() == 0) return;

    Node* cur = root;
    std::string remainKey = key;
    
    while(remainKey.size()>0)
    {
        char letter = remainKey[0];
        remainKey = remainKey.substr(1); //shorten the key
        
        typename std::vector<Node*>::iterator it = cur->children.begin();
        for(; it != cur->children.end(); it++)
        {
            if((*it)->label == letter) //child root already exist
            {
                cur = (*it);
                break;
            }
        }
        
        if(it == cur->children.end()) //children root doesn't exist
        {
            Node* p = new Node;
            p->label = letter;
            (cur->children).push_back(p);
            cur = p;
        }
    
    }
    (cur->val).push_back(value);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    std::vector<ValueType> empty;
    if(exactMatchOnly) //exact match
    {
        Node* cur = root;
        std::string remainKey = key;
        
        //go through each key of the character to find exact match character
        while(remainKey.size()>0)
        {
            char letter = remainKey[0];
            remainKey = remainKey.substr(1); //shorten the key
            
            typename std::vector<Node*>::iterator it = cur->children.begin();
            for(; it != cur->children.end(); it++)
            {
                if((*it)->label == letter) //child root already exist
                {
                    cur = (*it);
                    break;
                }
            }
            if(it == cur->children.end()) return empty;
        }
        return cur->val;
    }    
    else //single mismatching character
    {
        std::vector<ValueType> result;
        result = combineVectors(result, find(key, true));
        
        Node* cur = root;
        typename std::vector<Node*>::iterator it = cur->children.begin();
        for(; it != cur->children.end(); it++)
        {
            if((*it)->label == key[0])
            {
                cur = (*it);
                break;
            }
        }
        if(it == cur->children.end()) return empty;
        
        for(int x = 1; x < key.size(); x++)
        {
            //creating a single mismatch for each letter of the key
            std::string key_before = key.substr(0, x);
            std::string key_after = key.substr(x + 1);
            char key_this = key[x];
            
            Node* next = nullptr;
            typename std::vector<Node*>::iterator it2 = cur->children.begin();
            for(; it2 != cur->children.end(); it2++)
            {
                std::string newKey = key_before + (*it2)->label + key_after;
                //call recursion with one mismatch
                if((*it2)->label != key_this)
                    result = combineVectors(result, find(newKey, true));
                if((*it2)->label == key_this) next = (*it2);
            }
            
            //meaning next valid character doesnt exist, and two mismatches aren't allowed
            if(next == nullptr) return result;
            cur = next;
        }
        
        return result;
    }
    return empty;
}



//--------------------private functions-------------------
template <typename ValueType>
void Trie<ValueType>::clearUp(Node* cur){
    
    if(cur == nullptr) return;
    
    for(int x = 0; x < cur->children.size(); x++)
    {
        clearUp(cur->children[x]);
    }
    
    delete cur;
}

template <typename ValueType>
std::vector<ValueType> Trie<ValueType>::combineVectors(std::vector<ValueType> v1, std::vector<ValueType> v2) const
{
    std::vector<ValueType> result = v1;
    for(int x = 0; x < v2.size(); x++)
    {
        result.push_back(v2[x]);
    }
    
    return result;
}

#endif // TRIE_INCLUDED
