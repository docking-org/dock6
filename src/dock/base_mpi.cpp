#include <sstream>
#include <numeric>
#include "base_mpi.h"
#include "dockmol.h"
#include "trace.h"

using namespace std;

/************************************************/
// serializes a dockmol object into a string
// for passing with mpi
/************************************************/
void
dockmol_to_string(DOCKMol & mol, string & text)
{
    Trace trace( "::dockmol_to_string" );
    int             i;
    int num_residues_in_rec = mol.footprints.size(); 

    char tmp_char = '\n';

    ostringstream *ostr = new ostringstream();

    // write the molecule size information
    *ostr << mol.num_atoms << endl;
    *ostr << mol.num_bonds << endl;
    *ostr << mol.num_residues << endl;
    *ostr << num_residues_in_rec << endl;

    // Write the molecule text info
    ostr->write(mol.title.c_str(), mol.title.size());
    ostr->put(tmp_char);
    ostr->write(mol.mol_info_line.c_str(), mol.mol_info_line.size());
    ostr->put(tmp_char);
    ostr->write(mol.comment1.c_str(), mol.comment1.size());
    ostr->put(tmp_char);
    ostr->write(mol.comment2.c_str(), mol.comment2.size());
    ostr->put(tmp_char);
    ostr->write(mol.comment3.c_str(), mol.comment3.size());
    ostr->put(tmp_char);
    ostr->write(mol.energy.c_str(), mol.energy.size());
    ostr->put(tmp_char);
    ostr->write(mol.score_text_data.c_str(), mol.score_text_data.size());
    ostr->put(tmp_char);
    
    // write the atom information
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.x[i] << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.y[i] << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.z[i] << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.charges[i] << endl;
    for (i = 0; i < mol.num_atoms; i++) {
        ostr->write(mol.atom_types[i].c_str(), mol.atom_types[i].size());
        ostr->put(tmp_char);
    }
    for (i = 0; i < mol.num_atoms; i++) {
        ostr->write(mol.atom_names[i].c_str(), mol.atom_names[i].size());
        ostr->put(tmp_char);
    }
    for (i = 0; i < mol.num_atoms; i++) {
        ostr->write(mol.atom_residue_numbers[i].c_str(),
                    mol.atom_residue_numbers[i].size());
        ostr->put(tmp_char);
    }
    for (i = 0; i < mol.num_atoms; i++) {
        ostr->write(mol.subst_names[i].c_str(), mol.subst_names[i].size());
        ostr->put(tmp_char);
    }

    // hbond acceptor
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) 
        *ostr << mol.flag_acceptor[i] << endl;
    // hbond donator
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) 
        *ostr << mol.flag_donator[i] << endl;
    // hbond heavy atom acc 
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) 
        *ostr << mol.acc_heavy_atomid[i] << endl;

    // write the bond information
    for (i = 0; i < mol.num_bonds; i++)
        *ostr << mol.bonds_origin_atom[i] << endl;
    for (i = 0; i < mol.num_bonds; i++)
        *ostr << mol.bonds_target_atom[i] << endl;
    for (i = 0; i < mol.num_bonds; i++) {
        ostr->write(mol.bond_types[i].c_str(), mol.bond_types[i].size());
        ostr->put(tmp_char);
    }

    // write ring information
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.atom_ring_flags[i] << endl;
    for (i = 0; i < mol.num_bonds; i++)
        *ostr << mol.bond_ring_flags[i] << endl;

    // write the activation information
    *ostr << mol.num_active_atoms << endl;
    *ostr << mol.num_active_bonds << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.atom_active_flags[i] << endl;
    for (i = 0; i < mol.num_bonds; i++)
        *ostr << mol.bond_active_flags[i] << endl;

    // write molecular info YZ added to allow printing to header
    *ostr << mol.rot_bonds << endl;
    *ostr << mol.heavy_atoms << endl;
    *ostr << mol.mol_wt << endl;
    *ostr << mol.formal_charge << endl;

    // write the scoring information
    // modified by trent balius
    *ostr << mol.current_score << endl;
    stringstream ss;
    string sline;
    ss << mol.current_data;
    int count_current_data_lines = 0;

    while(getline(ss,sline,'\n')){
       count_current_data_lines++;
    }
    *ostr << count_current_data_lines << endl;

    ss.clear();

    ostr->write(mol.current_data.c_str(), mol.current_data.size());

    // write the footprints
    // added by trent balius
    for (i = 0; i < num_residues_in_rec;i++){ 
         *ostr << mol.footprints[i].resname << endl;
         *ostr << mol.footprints[i].resid << endl;
         *ostr << mol.footprints[i].vdw << endl;
         *ostr << mol.footprints[i].es << endl;
         *ostr << mol.footprints[i].hb << endl;
    }
    // hbond info is string a string on multiple lines.
    // added by trent balius
    sline = "";
    ss << mol.hbond_text_data;

    // count number of lines in hbond info
    int count_hbond_info_lines = 0;
    while(getline(ss,sline,'\n')){
       count_hbond_info_lines++;
    }
    *ostr << count_hbond_info_lines << endl;
    ostr->write(mol.hbond_text_data.c_str(), mol.hbond_text_data.size());
    ss.clear();

    // write solvation and color information
    for (i = 0; i < mol.num_atoms; i++) {
        ostr->write(mol.atom_color[i].c_str(), mol.atom_color[i].size());
        ostr->put(tmp_char);
    }
    *ostr << mol.total_dsol << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.atom_psol[i] << endl;
    for (i = 0; i < mol.num_atoms; i++)
        *ostr << mol.atom_apsol[i] << endl;

    // added by Guilherme D.
    #ifdef BUILD_DOCK_WITH_RDKIT
    //to pass through rdkit properties 

    ////to prepare the pns names
    std::vector<std::string> vecpnstmp{};
    vecpnstmp = mol.pns_name;
    std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string("")); 
    //// write RDKit info
    if (mol.fails_filt) {
       *ostr << 1 << endl;
    } else {
       *ostr << 0 << endl;
    }
    *ostr << mol.num_stereocenters << endl;
    *ostr << mol.num_spiro_atoms << endl;
    *ostr << mol.clogp << endl;
    *ostr << mol.tpsa << endl;
    *ostr << mol.num_arom_rings << endl;
    *ostr << mol.num_alip_rings << endl;
    *ostr << mol.num_sat_rings << endl;
    ostr->write(mol.smiles.c_str(), mol.smiles.size());
    ostr->put(tmp_char);
    *ostr << mol.qed_score << endl;
    *ostr << mol.sa_score << endl;
    *ostr << mol.esol << endl;
    *ostr << mol.pns << endl;
    if ( mol.pns == 0){
    }else{
        for (int z = 0; z < mol.pns; z++) {
            *ostr << mol.pns_name[z] << std::endl;
        }
    }
    *ostr << mol.MACCS_size << std::endl;
    for (int z = 0; z < (mol.MACCS_size-1); z++) {
        *ostr << mol.MACCS[z];
    }
    #endif

    *ostr << mol.bad_molecule << std::endl;
    // write amber information
    ostr->write(mol.amber_score_ligand_id.c_str(),
                mol.amber_score_ligand_id.size());
    ostr->put(tmp_char);
    trace.note( "amber_score_ligand_id = " + mol.amber_score_ligand_id );

    // process ostr into string
    text = ostr->str();
    delete          ostr;
    ostr = NULL;
}

/************************************************/
// repacks the string passed by mpi back to
// a dockmol object
/************************************************/
void
string_to_dockmol(DOCKMol & mol, string & text)
{
    Trace trace( "::string_to_dockmol" );
    //trace.note( "string = " + text );
    int             i;
    int num_residues_in_rec; 
    int             tmp_int;
    char            line[1000];

    mol.clear_molecule();

    istringstream  *istr = new istringstream(text.c_str());
   
    // read the molecule size information
    istr->getline(line, 1000);
    sscanf(line, "%d", &(mol.num_atoms));
    istr->getline(line, 1000);
    sscanf(line, "%d", &(mol.num_bonds));
    istr->getline(line, 1000);
    sscanf(line, "%d", &(mol.num_residues));

    // allocate the arrays
    mol.allocate_arrays(mol.num_atoms, mol.num_bonds, mol.num_residues);

    istr->getline(line, 1000);
    sscanf(line, "%d", &(num_residues_in_rec));

    // read the general molecule info
    istr->getline(line, 1000);
    mol.title = line;
    istr->getline(line, 1000);
    mol.mol_info_line = line;
    istr->getline(line, 1000);
    mol.comment1 = line;
    istr->getline(line, 1000);
    mol.comment2 = line;
    istr->getline(line, 1000);
    mol.comment3 = line;
    istr->getline(line, 1000);
    mol.energy = line;
    istr->getline(line, 1000);
    mol.score_text_data = line;

    // read the atom information
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.x[i]));
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.y[i]));
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.z[i]));
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.charges[i]));
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        mol.atom_types[i] = line;
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        mol.atom_names[i] = line;
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        mol.atom_residue_numbers[i] = line;
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        mol.subst_names[i] = line;
    }

    // hbond acceptor
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.flag_acceptor[i] = (tmp_int == 1) ? true : false;
    }
    // hbond donator
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.flag_donator[i] = (tmp_int == 1) ? true : false;
    }
    // hbond heavy atom acc
    // added by trent balius
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &(mol.acc_heavy_atomid[i]));
    }

    // read the bond information
    for (i = 0; i < mol.num_bonds; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &(mol.bonds_origin_atom[i]));
    }
    for (i = 0; i < mol.num_bonds; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &(mol.bonds_target_atom[i]));
    }
    for (i = 0; i < mol.num_bonds; i++) {
        istr->getline(line, 1000);
        mol.bond_types[i] = line;
    }

    // read ring information
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.atom_ring_flags[i] = (tmp_int == 1) ? true : false;
    }
    for (i = 0; i < mol.num_bonds; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.bond_ring_flags[i] = (tmp_int == 1) ? true : false;
    }

    // read the activation information
    istr->getline(line, 1000);
    sscanf(line, "%d", &(mol.num_active_atoms));
    istr->getline(line, 1000);
    sscanf(line, "%d", &(mol.num_active_bonds));
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.atom_active_flags[i] = (tmp_int == 1) ? true : false;
    }
    for (i = 0; i < mol.num_bonds; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%d", &tmp_int);
        mol.bond_active_flags[i] = (tmp_int == 1) ? true : false;
    }

    // read molecular info YZ added to allow printing to header
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.rot_bonds));
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.heavy_atoms));
    istr->getline(line, 1000);
    sscanf(line, "%f", &(mol.mol_wt));
    istr->getline(line, 1000);
    sscanf(line, "%f", &(mol.formal_charge));
   

    // read the scoring information
    istr->getline(line, 1000);
    sscanf(line, "%f", &(mol.current_score));

    // read number of current_data lines then the current_data
    istr->getline(line, 1000);
    int num_current_data = 0;
    sscanf(line, "%d",&(num_current_data));
    stringstream ss;
    for (i=0;i<num_current_data;i++){
        istr->getline(line, 1000);
        ss << line << endl;
    }
    mol.current_data = ss.str();
    ss.clear();
    trace.note( "current_data = " + mol.current_data );


    // read the footprints information
    // added by trent balius

    FOOTPRINT_ELEMENT temp_footprint_ele;
    mol.footprints.clear();

    for (i = 0; i < num_residues_in_rec;i++) {
        getline(*istr,temp_footprint_ele.resname);
        
        istr->getline(line, 1000);
        sscanf(line,"%d",&(temp_footprint_ele.resid));

        istr->getline(line, 1000);
        sscanf(line,"%lf",&(temp_footprint_ele.vdw));

        istr->getline(line, 1000);
        sscanf(line,"%lf",&(temp_footprint_ele.es));

        istr->getline(line, 1000);
        sscanf(line,"%d",&(temp_footprint_ele.hb));

        mol.footprints.push_back(temp_footprint_ele);
    }

    //get number of lines in hbond info
    // added by trent balius
    istr->getline(line, 1000);
    int num_lines_hbinfo;
    sscanf(line, "%d",&(num_lines_hbinfo));
    for (i=0;i<num_lines_hbinfo;i++){
        istr->getline(line, 1000);
        ss << line << endl;
    }
    mol.hbond_text_data = ss.str();
    ss.clear();

    // read color and solvation information
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        mol.atom_color[i] = line;
    }
    istr->getline(line, 1000);
    sscanf(line, "%f", &(mol.total_dsol));
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.atom_psol[i]));
    }
    for (i = 0; i < mol.num_atoms; i++) {
        istr->getline(line, 1000);
        sscanf(line, "%f", &(mol.atom_apsol[i]));
    }

    // Added by Guilherme D. (Sept/2020)
    #ifdef BUILD_DOCK_WITH_RDKIT
    //to pass through rdkit properties 
 
    ////to prepare the pains names 
    std::vector<std::string> vecpnstmp{};
    vecpnstmp = mol.pns_name;
    std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string("")); 

    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.fails_filt));  

    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.num_stereocenters));
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.num_spiro_atoms));
    istr->getline(line, 1000);
    sscanf(line, "%lf", &(mol.clogp));
    istr->getline(line, 1000);
    sscanf(line, "%lf", &(mol.tpsa));    
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.num_arom_rings)); 
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.num_alip_rings));     
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.num_sat_rings)); 
    istr->getline(line, 1000);
    mol.smiles = line;
    istr->getline(line, 1000);
    sscanf(line, "%lf", &(mol.qed_score));
    istr->getline(line, 1000);
    sscanf(line, "%lf", &(mol.sa_score));
    istr->getline(line, 1000);
    sscanf(line, "%lf", &(mol.esol));
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.pns));
    if (mol.pns == 0) {
    }else{
        for (int z = 1; z<=mol.pns; z++) {
            istr->getline(line, 1000);
            mol.pns_name.push_back(line); 
        }
    }
    //MACCS_size is required because when MACCS fp calc succeeded it outputs 167 chars
    //when it fails, it outputs 166 chars.
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.MACCS_size));
    boost::dynamic_bitset<> tmp_bitset;
    for (int z = 0; z<(mol.MACCS_size-1); z++) {
        istr->getline(line, 1000); 
        tmp_bitset.push_back(strtol(line,nullptr,2));
    }
    mol.MACCS = tmp_bitset;
    #endif
    //needed for mpi dbfiltering
    istr->getline(line, 1000);
    sscanf(line, "%u", &(mol.bad_molecule));


    // read amber information
    istr->getline(line, sizeof(line));
    mol.amber_score_ligand_id = line;
    trace.note( "amber_score_ligand_id = " + mol.amber_score_ligand_id ); 
    // free the array 
    delete          istr;
    istr = NULL;
}
/************************************************/
bool
Base_MPI::initialize_mpi(int *pargc, char ***pargv)
{

    // According to the MPI Forum: An MPI implementation is free to require
    // that the arguments in the C binding must be the arguments to main.
    // Consequently, this function expects pointers to argc and to argv.

    cout << "Initializing MPI Routines...\n";

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Init(pargc, pargv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    // For 1 processor mpi execution, turn off mpi and run as serial 
    if (comm_size == 1) {
        MPI_Finalize();
        return false;
    }
#endif

    mpi_work_unit = 1;
    active_clients = comm_size - 1;
    continue_mols = true;
    return true;
}

/************************************************/
void
Base_MPI::finalize_mpi()
{

#ifdef BUILD_DOCK_WITH_MPI
    cout << "Finalizing MPI Routines...\n";
    MPI_Finalize();
#endif

}

/************************************************/
bool
Base_MPI::continue_reading_mols()
{

    if (continue_mols == true) {
        return true;
    } else {
        continue_mols = true;
        return false;
    }

}

/************************************************/
bool
Base_MPI::is_master_node()
{
    return rank == 0;
}

/************************************************/
bool
Base_MPI::is_client_node()
{
    return rank != 0;
}

/************************************************/
void
Base_MPI::add_to_send_queue(DOCKMol & mol)
{
    DOCKMol         tmp_mol;

    tmp_mol.clear_molecule();
    send_queue.push_back(tmp_mol);
    copy_molecule(send_queue[send_queue.size() - 1], mol); 
    if (send_queue.size() >= mpi_work_unit) {
        continue_mols = false;
    }
}

/************************************************/
bool
Base_MPI::get_from_send_queue(DOCKMol & mol)
{
    if (recv_queue.size() > 0)
        recv_queue_request();
    if (send_queue.size() == 0)
        send_queue_request();
    mol.clear_molecule();
    if (send_queue.size() > 0) {

        copy_molecule(mol, send_queue[send_queue.size() - 1]);
        send_queue.pop_back();

        return true;
    } else
        return false;

}

/************************************************/
void
Base_MPI::add_to_recv_queue(RANKMol & mol)
{
    RANKMol         tmp_mol;

    tmp_mol.clear_molecule();
    recv_queue.push_back(tmp_mol);
    copy_molecule(recv_queue[recv_queue.size() - 1].mol, mol.mol);
    recv_queue[recv_queue.size() - 1].data = mol.data;
    recv_queue[recv_queue.size() - 1].score = mol.score;

}

/************************************************/
bool
Base_MPI::get_from_recv_queue(RANKMol & mol)
{
    if ((recv_queue.size() == 0) && (receive_flag == true)) {
        receive_recv_queue();
        receive_flag = false;
    }
    mol.clear_molecule();
    
    if (recv_queue.size() > 0) {
        copy_molecule(mol.mol, recv_queue[recv_queue.size() - 1].mol);
        mol.data = recv_queue[recv_queue.size() - 1].data;
        mol.score = recv_queue[recv_queue.size() - 1].score;

        recv_queue.pop_back();
        return true;
    }
    return false;

}

/************************************************/
bool
Base_MPI::listen_for_request()
{
#ifdef BUILD_DOCK_WITH_MPI
    MPI_Status      status;
#endif

    if (active_clients > 0) {
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(request, 2, MPI_INT, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD,
                 &status);
#endif
        return true;
    } else
        return false;

}

/************************************************/
bool
Base_MPI::is_receive_request()
{

    return request[1] == 2;

}

/************************************************/
bool
Base_MPI::is_send_request()
{

    return request[1] == 1;

}

/************************************************/
bool
Base_MPI::send_queue_request()
{
#ifdef BUILD_DOCK_WITH_MPI
    MPI_Status      status;
#endif
    DOCKMol         tmp_mol;
    int             i;
    char           *mols;
    string          str;
    int             size = 0;
    tmp_mol.clear_molecule();
    request[0] = rank;
    request[1] = 1;

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Send(request, 2, MPI_INT, 0, 100, MPI_COMM_WORLD);
    MPI_Recv(response, 2, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
#endif

    while (response[1] > 0) {

#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(&size, 1, MPI_INT, response[0], 100, MPI_COMM_WORLD, &status);
#endif
        mols = new char[size];
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(mols, size, MPI_UNSIGNED_CHAR, response[0], 100,
                 MPI_COMM_WORLD, &status);
#endif

        str = "";
        for (i = 0; i < size; i++)
            str += mols[i];

        send_queue.push_back(tmp_mol);
        string_to_dockmol(send_queue[send_queue.size() - 1], str);

        delete[]mols;
        mols = NULL;
        response[1]--;

    }

    if (send_queue.size() > 0)
        return true;
    else
        return false;

}

/************************************************/
void
Base_MPI::transmit_send_queue()
{
    response[0] = 0;
    response[1] = send_queue.size();

    if (response[1] == 0)
        active_clients--;

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Send(response, 2, MPI_INT, request[0], 100, MPI_COMM_WORLD);
#endif

    string str;
    int size;
    while (send_queue.size() > 0) {
        str = "";
        dockmol_to_string(send_queue[send_queue.size() - 1], str);
        send_queue.pop_back();

        size = str.size();
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Send(&size, 1, MPI_INT, request[0], 100, MPI_COMM_WORLD);
        MPI_Send((void *) str.c_str(), size, MPI_UNSIGNED_CHAR, request[0], 100,
                 MPI_COMM_WORLD);
#endif
    }

}

/************************************************/
void
Base_MPI::recv_queue_request()
{
#ifdef BUILD_DOCK_WITH_MPI
    MPI_Status      status;
#endif
    int             size;
    int             num_mols;
    string          str;
    float           score;
    string          data;

    request[0] = rank;
    request[1] = 2;

    response[0] = 0;
    response[1] = 0;

    num_mols = recv_queue.size();

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Send(request, 2, MPI_INT, 0, 100, MPI_COMM_WORLD);
    MPI_Recv(response, 2, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
    MPI_Send(&num_mols, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
    MPI_Send(stats, 3, MPI_INT, 0, 100, MPI_COMM_WORLD);
#endif

    while (recv_queue.size() > 0) {
        str = "";
        data = "";
        dockmol_to_string(recv_queue[recv_queue.size() - 1].mol, str);
        data = recv_queue[recv_queue.size() - 1].data;
        score = recv_queue[recv_queue.size() - 1].score;

        size = str.size();
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Send(&size, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
        MPI_Send((void *) str.c_str(), size, MPI_UNSIGNED_CHAR, 0, 100,
                 MPI_COMM_WORLD);
#endif

        size = data.size();
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Send(&size, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
        MPI_Send((void *) data.c_str(), size, MPI_UNSIGNED_CHAR, 0, 100,
                 MPI_COMM_WORLD);

        MPI_Send(&score, 1, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
#endif

        recv_queue.pop_back();
    }

}

/************************************************/
void
Base_MPI::receive_recv_queue()
{
#ifdef BUILD_DOCK_WITH_MPI
    MPI_Status      status;
#endif
    char           *text;
    RANKMol         tmp_mol;
    string          str;
    float           score = -MIN_FLOAT;  // arbitrarily large score
    string          data;
    int             i;
    int             num_mols = 0;
    int             size = 0;

    response[0] = 0;
    response[1] = 0;

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Send(response, 2, MPI_INT, request[0], 100, MPI_COMM_WORLD);
    MPI_Recv(&num_mols, 1, MPI_INT, request[0], 100, MPI_COMM_WORLD, &status);
    MPI_Recv(stats, 3, MPI_INT, request[0], 100, MPI_COMM_WORLD, &status);
#endif

    while (num_mols > 0) {
        recv_queue.push_back(tmp_mol);

#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(&size, 1, MPI_INT, request[0], 100, MPI_COMM_WORLD, &status);
#endif
        text = new char[size];
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(text, size, MPI_UNSIGNED_CHAR, request[0], 100, MPI_COMM_WORLD,
                 &status);
#endif
        str = "";
        for (i = 0; i < size; i++)
            str += text[i];
        delete[]text;
        text = NULL;
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(&size, 1, MPI_INT, request[0], 100, MPI_COMM_WORLD, &status);
#endif
        text = new char[size];

#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(text, size, MPI_UNSIGNED_CHAR, request[0], 100, MPI_COMM_WORLD,
                 &status);
#endif
        data = "";
        for (i = 0; i < size; i++)
            data += text[i];
        delete[]text;
        text = NULL;
#ifdef BUILD_DOCK_WITH_MPI
        MPI_Recv(&score, 1, MPI_FLOAT, request[0], 100, MPI_COMM_WORLD,
                 &status);
#endif

        string_to_dockmol(recv_queue[recv_queue.size() - 1].mol, str);
        recv_queue[recv_queue.size() - 1].data = data;
        recv_queue[recv_queue.size() - 1].score = score;

        num_mols--;
    }

}

/************************************************/
void
Base_MPI::barrier()
{

#ifdef BUILD_DOCK_WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

}
