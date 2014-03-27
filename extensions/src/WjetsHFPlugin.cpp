#include <WjetsHFPlugin.hpp>

#include <Processor.hpp>

#include <array>
#include <limits>


using namespace std;


WjetsHFPlugin::WjetsHFPlugin(string const &name) noexcept:
    Plugin(name)
{}


Plugin *WjetsHFPlugin::Clone() const
{
    return new WjetsHFPlugin(name);
}


void WjetsHFPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool WjetsHFPlugin::ProcessEvent()
{
    // Count different types of flavours
    unsigned n_b = 0, n_bbar = 0, n_c = 0, n_cbar = 0, n_ue = 0;
    
    for (ShowerParton const &p: (*reader)->GetShowerPartons())
    {
        int const absPdgId = abs(p.GetPdgId());
                        
        if (absPdgId != 4 and absPdgId != 5)  // skip all but b and c quarks
            continue;
        
        // Immediate daughters of beam particles       
        if (p.GetOrigin() == ShowerParton::Origin::Proton)
        {
            ++n_ue;
            continue;
        }
                
        if (p.GetPdgId() == 5) ++n_b;
        if (p.GetPdgId() == -5) ++n_bbar; //maybe 5 corresponds to bbar, doesn't matter
        
        if (p.GetPdgId() == 4) ++n_c;
        if (p.GetPdgId() == -4) ++n_cbar;        
    }
    
    //cout << n_b << ' ' << n_bbar << ' ' << n_c << ' '<< n_cbar << ' ' << n_ue << ' ';
    if (n_b != 0) //if there is a b quark in the event (=pair b)
    {
        if (n_b != n_bbar) throw runtime_error("WjetsHFPlugin::ProcessEvent: odd number of b-quarks.");
        else type = Type::W_qq;
    } 
    else if ((n_c == n_cbar) and (n_c != 0)) //if there is c and cbar (this is a pair!) in the shower
        type = Type::W_qq;
    
    else if ((n_c != 0) and (n_cbar != 0)) //if the number of c != the number of cbar: need to consider of there is a pair ccbar in the ME
    {
        unsigned n_c_me_out = 0, n_cbar_me_out = 0;
        int i = 0; // particles counter
        for (GenParticle const &p: (*reader)->GetHardGenParticles())
        {
            ++i;
            if (i>4) //for outgoing particles
            {
                if (p.GetPdgId() == 4) ++n_c_me_out;
                if (p.GetPdgId() == -4) ++n_cbar_me_out;
            }
        }
            
        if ((n_c_me_out > 0) and (n_cbar_me_out > 0)) // there are c and cbar in the ME final state: it is a pair
            type = Type::W_qq;
        else //no pair
            type = Type::W_c;
    }
    
    else if ((n_c != 0) or (n_cbar != 0)) // if only 1 type of c quark exists
        type = Type::W_c;
    
    else if (n_ue != 0) //heavy flavours which are immediate proton's daughters
        type = Type::W_other;
    
    else //no heavy flavours at all
        type = Type::W_light;
    
    return true;
}


WjetsHFPlugin::Type WjetsHFPlugin::GetDecision() const noexcept
{
    return type;
}