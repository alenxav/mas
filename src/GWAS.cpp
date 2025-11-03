// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
SEXP GWAS(Eigen::MatrixXf Y, Eigen::MatrixXf GEN, Eigen::MatrixXf M,
          int NumPCs = 3, int maxit = 500, float logtol=-8,
          int cores = 1, bool verb = true){
  
  // 
  // Null model - Multivariate SNP-BLUP
  // 
  
  // Start setup
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = GEN.cols(), f = M.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  if(verb){ Rcpp::Rcout << "Running absobtion of phenotypes (SY)\n"; }
  Eigen::MatrixXf y(n0,k), WM(n0,f), MU(f,k), iMM(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WM.col(j)=M.col(j).array()*Z.col(i).array();}
    iMM = (WM.transpose()*WM).inverse();
    MU.col(i) = iMM * WM.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WM*MU.col(i) ).array()*Z.col(i).array();
  }
  
  // Compute SX
  if(verb){  Rcpp::Rcout << "Running absobtion of marker scores (SZ)\n"; }
  iMM = (M.transpose()*M).inverse();
  for(int j=0; j<m; j++){ GEN.col(j) = (GEN.col(j) - M*(iMM*M.transpose()*GEN.col(j))).array(); }
  
  // Single value decomposition
  if(verb){ Rcpp::Rcout << "SVD of marker score matrix\n"; }
  Eigen::BDCSVD<Eigen::MatrixXf> svd(GEN, Eigen::ComputeThinU | Eigen::ComputeThinV );
  int npc;
  if(NumPCs>0){npc = NumPCs;}else{npc = svd.matrixU().cols();}
  Eigen::MatrixXf V = svd.matrixV().leftCols(npc);
  Eigen::VectorXf D = svd.singularValues().head(npc);
  Eigen::MatrixXf X = (svd.matrixU().leftCols(npc) * D.array().matrix().asDiagonal());
  int p = X.cols();
  
  // Sum of squares of X
  if(verb){ Rcpp::Rcout << "Set starting values for coefficients and variances\n"; }
  Eigen::MatrixXf XX(p,k); 
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  Eigen::VectorXf MSX = XX.colwise().sum().array(); MSX=MSX.array()*iN.array();
  // Variances
  iN = (n.array()-f).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSX.array()).matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  Eigen::VectorXf TrDinvXX(k);
  Eigen::MatrixXf Dinv(p,k);
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  // Bending
  float deflate = 1.0, deflateMax = 0.9;
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0;
  // Convergence control
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0;
  int numit = 0;
  // Prior for stability
  float df0 = 1.1;
  Eigen::MatrixXf Sb = vb*df0;
  Eigen::VectorXf Se = ve*df0;
  Eigen::VectorXf iNp = (n.array()+df0-f).inverse();
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Loop
  if(verb){ Rcpp::Rcout << "Starting Gauss-Seidel loop\n"; }
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    //for(int J=0; J<p; J++){
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVe.array();
      if(k==1){ b1 = RHS.array()/LHS.array(); }else{
        b1 = LHS.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(RHS);}
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    //ve = ve.array() / (n.array()-f);
    ve = (ve.array()+Se.array()) * iNp.array();
    iVe = ve.array().inverse();
    h2 = 1 - ve.array()/vy.array();
    
    // Get tilde-hat
    for(int i=0; i<k; i++){
      Dinv.col(i) = (XX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
      TrDinvXX(i)  = (XX.col(i).transpose() * Dinv.col(i));}
    TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXX(i)+df0);
        }else{ // Covariances
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrDinvXX(i)+TrDinvXX(j));
        }}}
    
    // Bending
    A = vb*1.0;
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= deflate;}}}
    if(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){  deflate -= 0.005;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}}
    EVDofA.compute(A);
    MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ 
      inflate = inflate - MinDVb*1.001;
      A.diagonal().array() += inflate;}
    iG = A.inverse();
    
    // Covariances
    cnv = log10((beta0.array()-b.array()).square().sum());  
    CNV1(numit) = cnv;
    CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
    CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
    if( std::isnan(cnv) && verb ){Rcpp::Rcout << "Numerical issue!!\n"; break;}
    
    // Print
    ++numit;
    if( numit % 100 == 0 && verb){ Rcpp::Rcout << "Iter: "<< numit << " || log10 Conv: "<< cnv << "\n"; } 
    if( numit > 10 && cnv<logtol && verb ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
    if( numit == maxit && verb ){ Rcpp::Rcout << "Model did not converge\n"; }
  }
  
  if(verb){ Rcpp::Rcout << "Fitting final model\n"; }
  // Fitting the model
  Eigen::MatrixXf hat = X*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( M * MU.col(i) + hat.col(i)).array(); }
  Eigen::MatrixXf beta = V*b;
  
  // Correlations
  if(verb){ Rcpp::Rcout << "Estimating correlations\n"; }
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  if(verb){  Rcpp::Rcout << "Convergence statistics\n"; }
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){ CNV1b(i)=CNV1(i);CNV2b(i)=CNV2(i);CNV3b(i)=CNV3(i);}
  
  // Null model Output
  Rcpp::List NullModelOutput = Rcpp::List::create(Rcpp::Named("FxdEffCoef")=MU,
                                                  Rcpp::Named("MarkerEffects")=beta,
                                                  Rcpp::Named("NumOfPCs")=p,
                                                  Rcpp::Named("FittedValues")=hat,
                                                  Rcpp::Named("Heritability")=h2,
                                                  Rcpp::Named("GenCorrelations")=GC,
                                                  Rcpp::Named("VarBeta")=vb,
                                                  Rcpp::Named("VarResiduals")=ve,
                                                  Rcpp::Named("ConvergenceBeta")=CNV1b,
                                                  Rcpp::Named("ConvergenceH2")=CNV2b,
                                                  Rcpp::Named("ConvergenceVar")=CNV3b,
                                                  Rcpp::Named("NumOfIterations")=numit);
  
  // 
  // Membership GWAS
  // 
  
  if(verb){ Rcpp::Rcout << "Starting Membership GWAS\n"; }
  
  // Membership matrix and other things
  int df;
  float tmp, VarE, dff, MSB;
  Eigen::MatrixXf ZM(n0,f), PZM(n0,f), ZPZ(f,f);
  Eigen::MatrixXf QtlVar(m,k), QtlHer(m,k), Gamma(m,f), Wald(m,k), PopWald(m,f);
  Eigen::VectorXf mum(f), ZPy(f), gamma(f), iC(n0), tmpCoef(n0);
  Eigen::VectorXf yZg(n0), PyZg(n0);
  Eigen::MatrixXf MSEx(m,k);
  
  // Loop across traits
  for(int trait=0; trait<k; trait++){
    
    // Don't run if the base model broke
    if( std::isnan(cnv) && verb ){Rcpp::Rcout << "GWAS CANNOT RUN!!\n"; break;}
    
    // degrees of freedom
    df = n(trait)-1.0; 
    dff = 1.0/(df-f);
    
    // Bypassing P matrix
    iC = 1.0 / ( D.array().square() + (ve(trait)/vb(trait,trait)) ).array();
    
    // Genome screening
    if(verb){ Rcpp::Rcout << "Genome-wide screening (Trait "<< (trait+1) << ")\n"; }
    for(int J=0; J<m; J++){
      
      // Print marker  
      if( J>0 && J % 1000 == 0 && verb ){ Rcpp::Rcout << "Screening marker: "<< J << "\n"; } 
      
      // Create Shizhong's Z matrix
      for(int i=0; i<f; i++){
        ZM.col(i) =  M.col(i).array() * GEN.col(J).array();
        tmpCoef = iC.array().matrix().asDiagonal() * X.transpose() * ZM.col(i);
        PZM.col(i) = ZM.col(i).array() - (X*tmpCoef).array();}
      
      // Compute coefficient
      ZPZ = ZM.transpose() * PZM; 
      ZPZ.diagonal().array() += 0.1; // Added for algorithmic stability
      ZPy = ZM.transpose() * e.col(trait);
      gamma = ZPZ.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(ZPy);
      Gamma.row(J) = gamma.array();
      
      // Variance components and heritability
      yZg = (y.col(trait).array()-(ZM*gamma).array())* Z.col(trait).array();
      PyZg = e.col(trait).array()-(PZM*gamma).array();
      VarE = (yZg.transpose()*PyZg);
      VarE = VarE*dff;
      tmp = (gamma.transpose() * ZPZ * gamma).trace();
      MSB = tmp/f;
      
      // Store main var components
      MSEx(J,trait) = VarE;
      if(tmp>0){
        Wald(J,trait) = tmp/VarE;
        if((MSB-VarE)>0){
          QtlVar(J,trait) = f*(MSB-VarE)/df;
          QtlHer(J,trait) = QtlVar(J,trait)/(QtlVar(J,trait)+VarE);
        }else{
          QtlVar(J,trait)=0.0;
          QtlHer(J,trait)=0.0;
        }
      }else{
        QtlVar(J,trait)=0.0;
        QtlHer(J,trait)=0.0;
        Wald(J,trait)=0.0;
      }
      PopWald.row(J) = gamma.array().square() * ZPZ.diagonal().array()/VarE;}   
    
    // End GWAS loop
  }
  
  // GWAS output
  Rcpp::List GwasModelOutput = Rcpp::List::create(Rcpp::Named("Wald")=Wald,
                                                  Rcpp::Named("WaldByPop")=PopWald,
                                                  Rcpp::Named("QTLvar")=QtlVar,
                                                  Rcpp::Named("QTLh2")=QtlHer,
                                                  Rcpp::Named("Effect")=Gamma,
                                                  Rcpp::Named("Membership")=M,
                                                  Rcpp::Named("MSE")=MSEx,
                                                  Rcpp::Named("N")=n);
  // Final output
  Rcpp::List OutputList = Rcpp::List::create(Rcpp::Named("POLY")=NullModelOutput, 
                                             Rcpp::Named("GWAS")=GwasModelOutput);
  OutputList.attr("class") = "GenModel";
  return OutputList;
  
}

// [[Rcpp::export]]
SEXP MLM(Eigen::MatrixXf Y,         // n x k responses (traits in columns)
         Eigen::MatrixXf X,         // n x f fixed-effect design
         Rcpp::List Z_list,         // list of R random-effect designs (n x p_r)
         int   maxit = 500,
         float logtol = -8.0,       // convergence threshold in log10 space
         int   cores = 1,
         bool  verbose = false,     // turn Rout printing on/off (default: off)
         float df0 = 1.1,           // prior df (MLM-style)
         bool  NonNegativeCorr = false,
         bool  InnerGS = false,     // single-sweep inner GS to solve k×k system
         bool  NoInv = false,       // avoid iG in marker update (use vb-based system)
         bool  XFA = false,         // apply Extended Factor Analytic reconstruction
         int   NumXFA = 3)          // number of top eigenpairs to keep in XFA
{
  if (cores != 1) Eigen::setNbThreads(cores);
  
  // Dimensions
  const int n = Y.rows();
  const int k = Y.cols();
  const int f = X.cols();
  const int R = Z_list.size();
  
  if (verbose) Rcpp::Rcout << "Rout: Start PEGSX | n=" << n << " k=" << k << " f=" << f << " R=" << R
                           << " | InnerGS=" << InnerGS
                           << " | NoInv=" << NoInv
                           << " | XFA=" << XFA
                           << " | NumXFA=" << NumXFA
                           << " | NonNegativeCorr=" << NonNegativeCorr << "\n";
  
  // Build mask W and zero-out NA in Y
  Eigen::MatrixXf W(n, k);
  for (int i = 0; i < n; ++i) {
    for (int t = 0; t < k; ++t) {
      if (std::isnan(Y(i, t))) {
        W(i, t) = 0.0f;
        Y(i, t) = 0.0f;
      } else {
        W(i, t) = 1.0f;
      }
    }
  }
  Eigen::VectorXf n_each = W.colwise().sum();          // counts per trait
  if (verbose) Rcpp::Rcout << "Rout: Built mask W and counts per trait\n";
  
  // Masked X per trait (WX) and iXX = (WX'WX)^+; estimate fixed effects b once (MLM style)
  std::vector<Eigen::MatrixXf> WX_list(k);
  std::vector<Eigen::MatrixXf> iXX_list(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(f, k);      // fixed-effect coefficients
  
  for (int t = 0; t < k; ++t) {
    WX_list[t].resize(n, f);
    for (int j = 0; j < f; ++j) WX_list[t].col(j) = X.col(j).array() * W.col(t).array();
    
    Eigen::MatrixXf XX = WX_list[t].transpose() * WX_list[t];
    iXX_list[t] = XX.completeOrthogonalDecomposition().pseudoInverse();
    
    Eigen::VectorXf RHS = WX_list[t].transpose() * Y.col(t);
    b.col(t) = iXX_list[t] * RHS;
  }
  if (verbose) Rcpp::Rcout << "Rout: Estimated fixed effects b (masked OLS per trait)\n";
  
  // Residuals after absorbing fixed effects: y = (Y - WX*b) masked by W
  Eigen::MatrixXf y(n, k);
  for (int t = 0; t < k; ++t)
    y.col(t) = (Y.col(t) - WX_list[t] * b.col(t)).array() * W.col(t).array();
  
  if (verbose) Rcpp::Rcout << "Rout: Computed masked residuals y after fixed-effect absorption\n";
  
  // Orthogonalize each Z_r to X (MLM-style): Zr <- Zr - X * ( (X'X)^+ X' Zr )
  Eigen::MatrixXf iXX_orth = (X.transpose() * X).completeOrthogonalDecomposition().pseudoInverse();
  std::vector<Eigen::MatrixXf> Zmats(R);
  std::vector<int> pR(R, 0);
  for (int r = 0; r < R; ++r) {
    Zmats[r] = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    pR[r] = Zmats[r].cols();
    if (verbose) Rcpp::Rcout << "Rout: Orthogonalizing Z[" << r << "] to X | p_r=" << pR[r] << "\n";
    
    for (int j = 0; j < pR[r]; ++j) {
      Eigen::VectorXf proj = iXX_orth * (X.transpose() * Zmats[r].col(j));
      Zmats[r].col(j) = Zmats[r].col(j) - X * proj;
    }
  }
  
  // Precompute per-effect ZZ(j,t) = sum_i Z(i,j)^2 * W(i,t), TrZSZ=sum over markers, and tilde=Z' y
  std::vector<Eigen::MatrixXf> ZZ_list(R);       // p_r x k
  std::vector<Eigen::VectorXf> TrZSZ_list(R);    // k
  std::vector<Eigen::MatrixXf> tilde_list(R);    // p_r x k
  
  for (int r = 0; r < R; ++r) {
    ZZ_list[r].resize(pR[r], k);
    for (int j = 0; j < pR[r]; ++j) {
      ZZ_list[r].row(j) = Zmats[r].col(j).array().square().matrix().transpose() * W;
    }
    TrZSZ_list[r] = ZZ_list[r].colwise().sum();
    tilde_list[r] = Zmats[r].transpose() * y;
    if (verbose) Rcpp::Rcout << "Rout: Precomputed ZZ and tilde for Z[" << r << "]\n";
  }
  
  // Initialize random coefficients u_r
  std::vector<Eigen::MatrixXf> u_list(R);
  for (int r = 0; r < R; ++r) u_list[r] = Eigen::MatrixXf::Zero(pR[r], k);
  
  // Residual variances ve initialization (with df0 prior used in subsequent updates)
  Eigen::VectorXf ssy   = y.colwise().squaredNorm();
  Eigen::VectorXf denom = (n_each.array() - f).matrix();
  for (int t = 0; t < k; ++t) denom(t) = std::max(denom(t), 1.0f);
  Eigen::VectorXf iN_mlm = denom.array().inverse().matrix();
  
  Eigen::VectorXf ve = (ssy.array() * iN_mlm.array()).matrix() * 0.5f;
  ve = ve.array().max(1e-8f);
  Eigen::VectorXf iVe = ve.array().inverse();
  if (verbose) Rcpp::Rcout << "Rout: Initialized ve from y sums-of-squares\n";
  
  // vb and iG per effect (MLM-like init): vb_diag(t) = ve(t) / (TrZSZ(t) * 1/(n-f))
  std::vector<Eigen::MatrixXf> vb_list(R);
  std::vector<Eigen::MatrixXf> iG_list(R);
  std::vector<float> bend_inflate(R, 0.0f);
  for (int r = 0; r < R; ++r) {
    vb_list[r].resize(k, k);
    vb_list[r].setZero();
    for (int t = 0; t < k; ++t) {
      float denom_r = std::max(TrZSZ_list[r](t) * iN_mlm(t), 1e-8f);
      vb_list[r](t, t) = ve(t) / denom_r;
    }
    iG_list[r] = vb_list[r].completeOrthogonalDecomposition().pseudoInverse();
    if (verbose) Rcpp::Rcout << "Rout: Initialized vb and iG for Z[" << r << "]\n";
  }
  
  // Priors: Sb (per-effect), Se, iNp
  std::vector<Eigen::MatrixXf> Sb_list(R);
  for (int r = 0; r < R; ++r) Sb_list[r] = vb_list[r] * df0;
  Eigen::VectorXf Se  = ve * df0;
  Eigen::VectorXf iNp = (n_each.array() + df0 - f).matrix();
  for (int t = 0; t < k; ++t) iNp(t) = 1.0f / std::max(iNp(t), 1.0f);
  if (verbose) Rcpp::Rcout << "Rout: Set priors Sb_list, Se, iNp (df0=" << df0 << ")\n";
  
  // Working residuals (start with y)
  Eigen::MatrixXf e = y;
  
  // RGS index per effect
  std::vector<std::vector<int>> RGS_index(R);
  for (int r = 0; r < R; ++r) {
    RGS_index[r].resize(pR[r]);
    for (int j = 0; j < pR[r]; ++j) RGS_index[r][j] = j;
  }
  
  // Convergence
  std::vector<Eigen::MatrixXf> u0_list(R);
  float cnv = 10.0f;
  int   numit = 0;
  std::random_device rd_dev;
  std::mt19937 gen(rd_dev());
  
  if (verbose) Rcpp::Rcout << "Rout: Starting GS iterations\n";
  
  // Iteration loop
  while (numit < maxit) {
    // Store previous u
    for (int r = 0; r < R; ++r) u0_list[r] = u_list[r];
    
    // Randomized Gauss–Seidel loops across effects
    for (int r = 0; r < R; ++r) {
      std::shuffle(RGS_index[r].begin(), RGS_index[r].end(), gen);
      
      const Eigen::MatrixXf& Zr  = Zmats[r];
      Eigen::MatrixXf&       ur  = u_list[r];
      Eigen::MatrixXf&       iGr = iG_list[r];
      const Eigen::MatrixXf& ZZr = ZZ_list[r];
      Eigen::MatrixXf&       vbr = vb_list[r];
      
      if (verbose) Rcpp::Rcout << "Rout: GS effect r=" << r << " begin | NoInv=" << NoInv << " InnerGS=" << InnerGS << "\n";
      
      for (int jj = 0; jj < pR[r]; ++jj) {
        int J = RGS_index[r][jj];
        
        Eigen::VectorXf u0 = ur.row(J);
        
        // diag terms from ZZ(J,:) .* iVe
        Eigen::VectorXf diagVec = (ZZr.row(J).transpose().array() * iVe.array()).matrix();
        
        Eigen::MatrixXf LHS(k, k);
        Eigen::VectorXf RHS(k);
        
        if (NoInv) {
          // LHS = vb * diag(diagVec); LHS(i,i) += 1
          LHS = vbr * diagVec.asDiagonal();
          for (int i = 0; i < k; ++i) LHS(i, i) += 1.0f;
          // RHS = vb * ( iVe .* (Z' e + ZZ .* u0) )
          Eigen::VectorXf base = (Zr.col(J).transpose() * e).array();
          base.array() += ZZr.row(J).transpose().array() * u0.array();
          base.array() *= iVe.array();
          RHS = vbr * base;
        } else {
          // LHS = iG + diag(diagVec)
          LHS = iGr;
          LHS.diagonal().array() += diagVec.array();
          // RHS = iVe .* (Z' e + ZZ .* u0)
          RHS = (Zr.col(J).transpose() * e).array();
          RHS.array() += ZZr.row(J).transpose().array() * u0.array();
          RHS.array() *= iVe.array();
        }
        
        // Solve for u1
        Eigen::VectorXf u1(k);
        if (InnerGS) {
          // Single sweep inner GS to solve LHS*u1 = RHS
          u1 = ur.row(J).transpose(); // start from current
          for (int i = 0; i < k; ++i) {
            float diag_ii = LHS(i, i);
            if (std::abs(diag_ii) < 1e-12f) continue;
            float accum = 0.0f;
            for (int l = 0; l < k; ++l) {
              if (l == i) continue;
              accum += LHS(i, l) * u1(l);
            }
            u1(i) = (RHS(i) - accum) / diag_ii;
          }
        } else {
          u1 = LHS.llt().solve(RHS);
        }
        
        ur.row(J) = u1.transpose();
        
        // Update residuals with mask
        e = (e - (Zr.col(J) * (u1 - u0).transpose()).cwiseProduct(W)).eval();
      }
      
      if (verbose) Rcpp::Rcout << "Rout: GS effect r=" << r << " end\n";
    }
    
    // Update residual variance ve with prior: ve = (sum(e*y) + Se) * iNp
    Eigen::VectorXf new_ve = (e.cwiseProduct(y)).colwise().sum();
    new_ve = (new_ve.array() + Se.array()) * iNp.array();
    ve  = new_ve.array().max(1e-8f);
    iVe = ve.array().inverse();
    
    if (verbose) Rcpp::Rcout << "Rout: Updated residual variances ve (with priors)\n";
    
    // Update vb_r (per effect) and iG_r; apply NonNegativeCorr and XFA if requested; then bend
    for (int r = 0; r < R; ++r) {
      // TildeHat = u_r' * (Z_r' y)
      Eigen::MatrixXf TildeHat = u_list[r].transpose() * tilde_list[r];
      Eigen::MatrixXf& vbr = vb_list[r];
      
      // Diagonal (with Sb prior and df0)
      for (int i = 0; i < k; ++i) {
        float denom_i = std::max(TrZSZ_list[r](i), 1e-8f);
        vbr(i,i) = (TildeHat(i,i) + Sb_list[r](i,i)) / (denom_i + df0);
      }
      // Off-diagonal (no prior)
      for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j) {
          float denom_ij = std::max(TrZSZ_list[r](i) + TrZSZ_list[r](j), 1e-8f);
          float cov_ij   = (TildeHat(i,j) + TildeHat(j,i)) / denom_ij;
          vbr(i,j) = cov_ij;
          vbr(j,i) = cov_ij;
        }
        
        // Enforce non-negative covariances if requested
        if (NonNegativeCorr) {
          vbr = vbr.array().cwiseMax(0.0f).matrix();
          vbr = 0.5f * (vbr + vbr.transpose());
          if (verbose) Rcpp::Rcout << "Rout: Applied NonNegativeCorr to vb for effect r=" << r << "\n";
        }
        
        // Optional XFA: reconstruct GC using top NumXFA eigenpairs, then scale back to vb
        if (XFA && NumXFA > 0) {
          // Compute GC from current vb
          Eigen::VectorXf sd = vbr.diagonal().array().sqrt();
          for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
          Eigen::VectorXf inv_sd = sd.array().inverse();
          Eigen::MatrixXf GC = inv_sd.asDiagonal() * vbr * inv_sd.asDiagonal();
          
          // EVD of GC
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> esGC(GC);
          Eigen::VectorXf evals = esGC.eigenvalues();
          Eigen::MatrixXf evecs = esGC.eigenvectors(); // columns are eigenvectors
          
          // Reconstruct using top NumXFA eigenpairs (largest eigenvalues)
          int useF = std::min(NumXFA, k);
          Eigen::MatrixXf GC_red = Eigen::MatrixXf::Zero(k, k);
          for (int ii = 0; ii < useF; ++ii) {
            int idx = k - 1 - ii; // largest
            float lam = evals(idx);
            Eigen::VectorXf v = evecs.col(idx);
            GC_red += lam * (v * v.transpose());
          }
          // Set diagonal to 1 (correlations)
          for (int i = 0; i < k; ++i) GC_red(i, i) = 1.0f;
          
          // Convert back to covariance using current diagonals
          for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
              vbr(i, j) = GC_red(i, j) * std::sqrt(std::max(vbr(i, i), 1e-12f) * std::max(vbr(j, j), 1e-12f));
            }
          }
          if (verbose) Rcpp::Rcout << "Rout: Applied XFA with NumXFA=" << useF << " for effect r=" << r << "\n";
        }
        
        // Bending if needed
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(vbr);
        float min_ev = es.eigenvalues().minCoeff();
        if (min_ev < 0.001f) {
          float need = std::abs(min_ev * 1.1f);
          bend_inflate[r] = std::max(bend_inflate[r], need);
        }
        vbr.diagonal().array() += bend_inflate[r];
        
        // Precision (computed regardless of NoInv)
        iG_list[r] = vbr.completeOrthogonalDecomposition().pseudoInverse();
        
        if (verbose) Rcpp::Rcout << "Rout: Updated vb and iG for Z[" << r << "], bend=" << bend_inflate[r] << "\n";
    }
    
    // Convergence on u changes (max across effects)
    double ss_max = 0.0;
    for (int r = 0; r < R; ++r) {
      double diff = (u0_list[r].array() - u_list[r].array()).square().sum();
      ss_max = std::max(ss_max, diff);
    }
    cnv = std::log10(static_cast<float>(std::max(ss_max, 1e-20)));
    
    ++numit;
    if (verbose) Rcpp::Rcout << "Rout: Iter " << numit << " | log10(cnv)=" << cnv << " | ss_max=" << ss_max << "\n";
    
    // Require at least a few iterations before stopping
    if (numit >= 5 && (cnv < logtol || std::isnan(cnv))) {
      if (verbose) Rcpp::Rcout << "Rout: Convergence reached or NaN detected, breaking loop\n";
      break;
    }
  }
  
  if (verbose) Rcpp::Rcout << "Rout: Iterations finished. numit=" << numit << " cnv=" << cnv << "\n";
  
  // Fitted values: hat = X*b + sum_r Z_r * u_r
  Eigen::MatrixXf hat = X * b;
  for (int r = 0; r < R; ++r) hat.noalias() += Zmats[r] * u_list[r];
  if (verbose) Rcpp::Rcout << "Rout: Computed fitted values hat\n";
  
  // Genetic correlations per effect (standardize vb)
  Rcpp::List GC_out(R);
  Rcpp::List vb_out(R);
  Rcpp::List u_out(R);
  Rcpp::NumericVector bend_out(R);
  
  for (int r = 0; r < R; ++r) {
    Eigen::VectorXf sd = vb_list[r].diagonal().array().sqrt();
    for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
    Eigen::VectorXf inv_sd = sd.array().inverse();
    Eigen::MatrixXf GC = inv_sd.asDiagonal() * vb_list[r] * inv_sd.asDiagonal();
    
    vb_out[r]   = vb_list[r];
    GC_out[r]   = GC;
    u_out[r]    = u_list[r];
    bend_out[r] = bend_inflate[r];
  }
  if (verbose) Rcpp::Rcout << "Rout: Prepared outputs vb_list, GC_list, u, bend\n";
  
  return Rcpp::List::create(
    Rcpp::Named("b")        = b,         // f x k fixed-effect coefficients
    Rcpp::Named("u")        = u_out,     // list of p_r x k random-effect coefficients
    Rcpp::Named("vb_list")  = vb_out,    // list of k x k covariances per random effect
    Rcpp::Named("GC_list")  = GC_out,    // list of k x k genetic correlations per random effect
    Rcpp::Named("ve")       = ve,        // k residual variances
    Rcpp::Named("hat")      = hat,       // n x k fitted values
    Rcpp::Named("its")      = numit,     // iterations
    Rcpp::Named("cnv")      = cnv,       // log10 convergence score (max across effects)
    Rcpp::Named("bend")     = bend_out   // per-effect bending magnitude
  );
}

