/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    paraboloidVelocityFvPatchVectorField

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    paraboloidVelocityFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/
/*
#ifndef…#define…#endif：这是经典的“包含保护（include guard）”，防止头文件被多次包含导致重定义错误。
*/
#ifndef paraboloidVelocityFvPatchVectorField_H
#define paraboloidVelocityFvPatchVectorField_H

#include "fvPatchFields.H"//定义了所有 fvPatchField 的基类、patch 网格结构等。
#include "fixedValueFvPatchFields.H"//提供了 “固定值” 边界条件的具体实现类 fixedValueFvPatchVectorField，我们的新类paraboloidVelocityFvPatchVectorField继承自它。

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam//namespace Foam：所有 OpenFOAM 代码都放在 Foam 命名空间里，防止与用户或第三方代码冲突
{

/*---------------------------------------------------------------------------*\
              Class paraboloidVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class paraboloidVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField//paraboloidVelocityFvPatchVectorField 继承自 fixedValueFvPatchVectorField
{
    // Private data

        //- Peak velocity magnitude
        
        scalar Umax_;
                
        scalar Rm_;
        
        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;
        vector z_;
        vector ctr_;
// 不写修饰符 默认是私有属性/ 这些都是本类特有的“状态”（state），用来在 updateCoeffs() 中计算抛物线流型的各个几何／速度参数

public:

    //- Runtime type information
    TypeName("paraboloidVelocity");//这是 OpenFOAM 自己的宏，注册这个类的名字，以便在字典里写： 边界条件 type            paraboloidVelocity;


    // Constructors

        //- Construct from patch and internal field
        paraboloidVelocityFvPatchVectorField//最基本的构造：给定 patch（边界面）和内场（例如 U、C 等），构造一个默认行为的边界场。
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        paraboloidVelocityFvPatchVectorField//这构造方法额外读取运行时字典（fvPatchFieldDict）里的参数，比如 Umax, Rm, n, ctr，并存入私有成员。
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given paraboloidVelocityFvPatchVectorField
        //  onto a new patch
        paraboloidVelocityFvPatchVectorField//用于在动态网格、重分区（decomposePar/reconstructPar）或 patch 重用时，把一个已存在的边界条件“映射（map）”到新的 patch 上。
        (
            const paraboloidVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const//返回一个 new 出来的深拷贝，配合工厂和多态使用，保证任何时候都能生成本类型的实例。
        {
            return tmp<fvPatchVectorField>
            (
                new paraboloidVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        paraboloidVelocityFvPatchVectorField
        (
            const paraboloidVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new paraboloidVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
               
        scalar& Umax()//
        {
            return Umax_;
        }
        
        scalar& Rm()
        {
            return Rm_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }
        
        vector& z()
        {
            return z_;
        }        
     
        vector& ctr()
        {
            return ctr_;
        }

        //- Update coefficients
        virtual void updateCoeffs();
        //这是 fvPatchField 系列标准的“更新系数”接口，OpenFOAM 每次求解迭代都会调用它。
        //在这里，用户要根据 Umax_、Rm_、n_、ctr_ 等私有数据，计算出每个 patch face 上的实际速度矢量，并通过 operator== 赋值。


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
