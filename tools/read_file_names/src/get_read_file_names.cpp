#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
bool isALeaf(IBank* bank)
{
    if(dynamic_cast<BankComposite*> (bank) == 0) return true;
    return false;
}

/********************************************************************************/
void printCis(IBank* bank, int* value_i, size_t depth=0)
{
    LOCAL (bank);
    if (depth==1){
        cout <<"C_"<<*value_i<<" ";
        (*value_i)+=1;
    }
    if (isALeaf(bank)){
        cout << bank->getId() << " ";
    }


    BankComposite* composite = dynamic_cast<BankComposite*> (bank);
    if (composite != 0)
    {
        const std::vector<IBank*> children = composite->getBanks();
        for (size_t i=0; i<children.size(); i++)  {  printCis (children[i], value_i, depth+1); }
    }

    if (depth==1) cout<<endl;
}

/********************************************************************************/
void dump (IBank* bank, size_t depth=0)
{
    LOCAL (bank);

    for (size_t i=0; i<depth; i++)  { cout << "  "; }  cout << bank->getId() << endl;

    BankComposite* composite = dynamic_cast<BankComposite*> (bank);
    if (composite != 0)
    {
        const std::vector<IBank*> children = composite->getBanks();
        for (size_t i=0; i<children.size(); i++)  {  dump (children[i], depth+1); }
    }
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankDump");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We dump the bank hierarchy. */
        //dump (Bank::open (options->getStr(STR_URI_INPUT)));

        /** We dump the bank leaves hierarchy. */
        int value_i=1;
        printCis (Bank::open (options->getStr(STR_URI_INPUT)), &value_i);
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
