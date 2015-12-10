﻿/*
 * Original author: Nicholas Shulman <nicksh .at. u.washington.edu>,
 *                  MacCoss Lab, Department of Genome Sciences, UW
 *
 * Copyright 2011 University of Washington - Seattle, WA
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
using System;
using System.Collections.Generic;
using System.Linq;

namespace pwiz.Common.DataBinding
{
    public interface IFilterOperation
    {
        string OpName { get; }
        string DisplayName { get; }
        bool IsValidFor(ColumnDescriptor columnDescriptor);
        Type GetOperandType(ColumnDescriptor columnDescriptor);
        Predicate<object> MakePredicate(ColumnDescriptor columnDescriptor, string operand);
    }

    public static class FilterOperations
    {

        public static readonly IFilterOperation OP_HAS_ANY_VALUE = new FilterOperation("", "Has Any Value", null,
                                                                         (cd, operand) => rowNode => true);

        public static readonly IFilterOperation OP_EQUALS 
            = new FilterOperation("equals", "Equals", typeof (object), FnEquals);

        public static readonly IFilterOperation OP_NOT_EQUALS 
            = new FilterOperation("<>", "Does Not Equal", typeof (object), FnNotEquals);

        public static readonly IFilterOperation OP_IS_BLANK
            = new FilterOperation("isnullorblank", "Is Blank", null, FnIsBlank);

        public static readonly IFilterOperation OP_IS_NOT_BLANK
            = new FilterOperation("isnotnullorblank", "Is Not Blank", null, FnIsNotBlank);

        public static readonly IFilterOperation OP_IS_GREATER_THAN
            = new FilterOperation(">", "Is Greater Than", typeof (object), MakeFnCompare(i => i > 0));
        public static readonly IFilterOperation OP_IS_LESS_THAN
            = new FilterOperation("<", "Is Less Than", typeof(object), MakeFnCompare(i => i < 0));
        public static readonly IFilterOperation OP_IS_GREATER_THAN_OR_EQUAL
            = new FilterOperation(">=", "Is Greater Than Or Equal To", typeof(object), MakeFnCompare(i => i >= 0));
        public static readonly IFilterOperation OP_IS_LESS_THAN_OR_EQUAL
            = new FilterOperation("<=", "Is Less Than Or Equal To", typeof(object), MakeFnCompare(i => i <= 0));

        public static readonly IFilterOperation OP_CONTAINS = new StringFilterOperation("contains", "Contains", FnContains);
        public static readonly IFilterOperation OP_NOT_CONTAINS = new StringFilterOperation("notcontains", "Does Not Contain", FnNotContains);
        public static readonly IFilterOperation OP_STARTS_WITH = new StringFilterOperation("startswith", "Starts With", FnStartsWith);
        public static readonly IFilterOperation OP_NOT_STARTS_WITH = new StringFilterOperation("notstartswith", "Does Not Start With", FnNotStartsWith);

        private static readonly IList<IFilterOperation> LST_FILTER_OPERATIONS = Array.AsReadOnly(new[]
                                                                                           {
                                                                                               OP_HAS_ANY_VALUE,
                                                                                               OP_EQUALS,
                                                                                               OP_NOT_EQUALS,
                                                                                               OP_IS_BLANK,
                                                                                               OP_IS_NOT_BLANK,
                                                                                               OP_IS_GREATER_THAN,
                                                                                               OP_IS_LESS_THAN,
                                                                                               OP_IS_GREATER_THAN_OR_EQUAL,
                                                                                               OP_IS_LESS_THAN_OR_EQUAL,
                                                                                               OP_CONTAINS,
                                                                                               OP_NOT_CONTAINS,
                                                                                               OP_STARTS_WITH,
                                                                                               OP_NOT_STARTS_WITH
                                                                                           });

        private static readonly IDictionary<string, IFilterOperation> DICT_FILTER_OPERATIONS =
            LST_FILTER_OPERATIONS.ToDictionary(op => op.OpName, op => op);
        public static IFilterOperation GetOperation(string name)
        {
            IFilterOperation result;
            DICT_FILTER_OPERATIONS.TryGetValue(name, out result);
            return result;
        }
        public static IList<IFilterOperation> ListOperations()
        {
            return LST_FILTER_OPERATIONS;
        }

        public static Predicate<object> FnEquals(ColumnDescriptor columnDescriptor, string strOperand)
        {
            object operand = ConvertOperand(columnDescriptor, strOperand);
            return value => Equals(ConvertValue(columnDescriptor, value), operand);
        }
        public static Predicate<object> FnNotEquals(ColumnDescriptor columnDescriptor, string strOperand)
        {
            object operand = ConvertOperand(columnDescriptor, strOperand);
            return value => !Equals(ConvertValue(columnDescriptor, value), operand);
        }
        public static Predicate<object> FnIsBlank(ColumnDescriptor columnDescriptor, string operand)
        {
            return value => null == value || Equals(value, "");
        }
        public static Predicate<object> FnIsNotBlank(ColumnDescriptor columnDescriptor, string operand)
        {
            return value => null != value && !Equals(value, "");
        }
        public static Func<ColumnDescriptor, string, Predicate<object>> MakeFnCompare(Predicate<int> comparisonPredicate)
        {
            return (columnDescriptor, strOperand)
                   =>
                       {
                           object operand = ConvertOperand(columnDescriptor, strOperand);
                           return value => null != value && comparisonPredicate(columnDescriptor.DataSchema.Compare(value, operand));
                       };
            
        }
        public static Predicate<object> FnContains(ColumnDescriptor columnDescriptor, string strOperand)
        {
            return value => null != value && value.ToString().IndexOf(strOperand, StringComparison.Ordinal) >= 0;
        }
        public static Predicate<object> FnNotContains(ColumnDescriptor columnDescriptor, string strOperand)
        {
            return value => null != value && value.ToString().IndexOf(strOperand, StringComparison.Ordinal) < 0;
        }
        public static Predicate<object> FnStartsWith(ColumnDescriptor columnDescriptor, string strOperand)
        {
            return value => null != value && value.ToString().StartsWith(strOperand);
        }
        public static Predicate<object> FnNotStartsWith(ColumnDescriptor columnDescriptor, string strOperand)
        {
            return value => null != value && !value.ToString().StartsWith(strOperand);
        }


        public static object ConvertOperand(ColumnDescriptor columnDescriptor, string operand)
        {
            var type = columnDescriptor.WrappedPropertyType;
            if (typeof(char) == type)
            {
                if (operand.Length != 1)
                {
                    return null;
                }
                return operand[0];
            }
            if (typeof(bool) == type)
            {
                return Convert.ToBoolean(operand);
            }
            if (type.IsEnum)
            {
                if (string.IsNullOrEmpty(operand))
                {
                    return null;
                }
                try
                {
                    return Enum.Parse(type, operand);
                }
                catch
                {
                    return Enum.Parse(type, operand, true);
                }
            }
            if (type.IsPrimitive)
            {
                if (typeof(double) == type)
                {
                    return Convert.ToDouble(operand);
                }
                if (typeof(float) == type)
                {
                    return Convert.ToSingle(operand);
                }
                if (typeof(decimal) == type)
                {
                    return Convert.ToDecimal(operand);
                }
                return Convert.ToInt64(operand);
            }
            return Convert.ToString(operand);
        }
        public static object ConvertValue(ColumnDescriptor columnDescriptor, object value)
        {
            var type = columnDescriptor.WrappedPropertyType;
            if (typeof(char) == type)
            {
                return value as char?;
            }
            if (typeof(bool) == type)
            {
                return value as bool?;
            }
            if (type.IsEnum)
            {
                return value;
            }
            if (type.IsPrimitive)
            {
                if (typeof(double) == type)
                {
                    return Convert.ToDouble(value);
                }
                if (typeof(float) == type)
                {
                    return Convert.ToSingle(value);
                }
                if (typeof(decimal) == type)
                {
                    return Convert.ToDecimal(value);
                }
                return Convert.ToInt64(value);
            }
            return Convert.ToString(value);
        }

        class FilterOperation : IFilterOperation
        {
            private readonly Func<ColumnDescriptor, string, Predicate<object>> _fnMakePredicate;
            private readonly Type _operandType;
            public FilterOperation(string opName, string displayName, Type operandType, Func<ColumnDescriptor, string, Predicate<object>> fnMakePredicate)
            {
                OpName = opName;
                DisplayName = displayName;
                _operandType = operandType;
                _fnMakePredicate = fnMakePredicate;
            }
            public string OpName { get; private set; }
            public string DisplayName { get; private set; }
            public virtual bool IsValidFor(ColumnDescriptor columnDescriptor)
            {
                return true;
            }
            public virtual Type GetOperandType(ColumnDescriptor columnDescriptor)
            {
                return _operandType;
            }
            public Predicate<object> MakePredicate(ColumnDescriptor columnDescriptor, string operand)
            {
                return _fnMakePredicate(columnDescriptor, operand);
            }
        }
        class StringFilterOperation : FilterOperation
        {
            public StringFilterOperation(string opName, string displayName, Func<ColumnDescriptor, string, Predicate<object>> fnMakePredicate) : base(opName, displayName, typeof(string), fnMakePredicate)
            {
            }
            public override bool IsValidFor(ColumnDescriptor columnDescriptor)
            {
                return typeof (string) == columnDescriptor.WrappedPropertyType;
            }
        }
    }
}